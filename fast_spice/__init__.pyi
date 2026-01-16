# fast_spice/__init__.pyi
from typing import List, Dict, Union, Any, Optional

class Circuit:
    """
    Classe principal do simulador Fast SPICE.
    Permite a construção e simulação de circuitos de eletrônica de potência com 
    suporte a Steady-State (Shooting Method) e Transiente.
    """

    # Atributos de estado expostos pelo motor C++
    has_state: bool 
    """Indica se o circuito possui um estado interno (V no Capacitor, I no Indutor) salvo."""
    
    last_state: List[float] 
    """Vetor de estado da última simulação realizada."""

    def __init__(self) -> None: 
        """Inicializa uma nova netlist vazia."""
        ...

    def set_params(self, max_iter: int, tol: float, eps: float = 1e-6) -> None:
        """
        Configura os parâmetros numéricos do Solver.

        Args:
            max_iter (int): Número máximo de iterações do Newton-Raphson para o Shooting Method.
            tol (float): Tolerância de convergência para o erro de regime permanente.
            eps (float): Perturbação infinitesimal para o cálculo numérico do Jacobiano (Sensibilidade).
        """
        ...

    def reset_state(self) -> None:
        """
        Limpa o estado interno salvo (last_state), forçando a próxima simulação 
        a começar do zero (Cold Start).
        """
        ...

    def update_value(self, name: str, val: float) -> None:
        """
        Atualiza o valor de um componente existente na netlist sem precisar reconstruí-la.
        Útil para varreduras de parâmetros ou otimização.

        Args:
            name: ID do componente (ex: "R1").
            val: Novo valor numérico.
        """
        ...

    # --- Métodos de Construção da Netlist ---

    def add_resistor(self, name: str, n1: str, n2: str, val: float) -> None:
        """Adiciona um resistor (Ohms) entre os nós n1 e n2."""
        ...

    def add_capacitor(self, name: str, n1: str, n2: str, val: float) -> None:
        """Adiciona um capacitor (Farads) entre os nós n1 e n2."""
        ...

    def add_inductor(self, name: str, n1: str, n2: str, val: float) -> None:
        """Adiciona um indutor (Henrys) entre os nós n1 e n2."""
        ...

    def add_dc_source(self, name: str, p: str, n: str, val: float) -> None:
        """Adiciona uma fonte de tensão DC (Volts) entre o polo positivo (p) e negativo (n)."""
        ...

    def add_mosfet(self, name: str, d: str, s: str) -> None:
        """
        Adiciona uma chave controlada ideal (MOSFET/IGBT) entre Dreno (d) e Source (s).
        Nota: Este método adiciona automaticamente um diodo antiparalelo à chave.
        """
        ...

    def add_diode(self, name: str, a: str, k: str) -> None:
        """Adiciona um diodo ideal entre Anodo (a) e Catodo (k)."""
        ...

    def add_transformer(self, name: str, p1: str, p2: str, s1: str, s2: str, ratio: float) -> None:
        """
        Adiciona um transformador ideal.
        ratio: Relação de espiras (N1/N2).
        """
        ...

    def set_pwm(
        self, 
        sw_name: str, 
        freq: float, 
        duty: float, 
        phase: float = 0.0, 
        dead_time: float = 0.0
    ) -> None:
        """
        Configura o sinal de gate para uma chave específica.

        Args:
            sw_name: Nome da chave definida em add_mosfet.
            freq: Frequência de chaveamento (Hz).
            duty: Razão cíclica (0.0 a 1.0).
            phase: Atraso de fase em relação ao período (0.0 a 1.0).
            dead_time: Tempo morto (segundos) aplicado ao início da condução.
        """
        ...

    def run(
        self, 
        mode: str, 
        dt: float, 
        t_end: float, 
        period: float = 0.0, 
        probes: List[str] = [],
        export_matrices: bool = False
    ) -> Dict[str, Any]:
        """
        Executa o motor de simulação.

        Args:
            mode: "transient" (simulação no tempo) ou "shooting" (busca de regime permanente).
            dt: Passo de tempo para integração (segundos).
            t_end: Tempo total para coleta de dados após convergência (segundos). Se 0, apenas o estado final é retornado.
            period: Período fundamental do circuito (obrigatório para modo "shooting").
            probes: Lista de sinais a monitorar, ex: ["V(out)", "I(L1)"].
            export_matrices: Se True, o dicionário de retorno incluirá o Jacobiano e sua inversa.

        Returns:
            Um dicionário contendo:
                - 't': Lista de tempos.
                - 'signals': Dict com as formas de onda dos probes.
                - 'stats': Estatísticas (mean, rms, max, min, ripple) de cada sinal.
                - 'jacobian': Matriz J do sistema (se export_matrices=True).
                - 'inverse_jacobian': Matriz B do sistema (se export_matrices=True).
        """
        ...