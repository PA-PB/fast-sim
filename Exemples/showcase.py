import numpy as np
import matplotlib.pyplot as plt
import fast_spice
import time

class LLCFullAnalysis:
    def __init__(self):
      
        self.params = {
            "Vin": 400.0,
            "Cr": 24e-9,      # Capacitor Ressonante (24nF)
            "Lr": 60e-6,      # Indutor Ressonante (60uH)
            "Lm": 280e-6,     # Indutância de Magnetização (280uH)
            "n_ratio": 8.2,   # Relação de Espiras
            "Co": 200e-6,     # Capacitor de Saída
            "Rload": 5.0,     # Carga
            "dead_time": 200e-9
        }
        
        # Frequência de Ressonância Teórica
        self.fr = 1 / (2 * np.pi * np.sqrt(self.params["Lr"] * self.params["Cr"]))
        print(f"Frequência de Ressonância: {self.fr/1e3:.2f} kHz")

        self.circuit = fast_spice.Circuit()
        self._build_netlist()

    def _build_netlist(self):
        p = self.params
        c = self.circuit
        
        # Topologia
        c.add_dc_source("VDC", "vdc", "gnd", p["Vin"])
        c.add_mosfet("QH", "vdc", "sw")
        c.add_mosfet("QL", "sw", "gnd")
        
        c.add_capacitor("Cr", "sw", "tank", p["Cr"])
        c.add_inductor("Lr", "tank", "pri", p["Lr"])
        c.add_inductor("Lm", "pri", "gnd", p["Lm"])
        
        c.add_transformer("Tx", "pri", "gnd", "sec1", "sec2", p["n_ratio"])
        c.add_diode("D1", "sec1", "out")
        c.add_diode("D2", "sec2", "out")
        c.add_diode("D3", "gnd", "sec1")
        c.add_diode("D4", "gnd", "sec2")
        
        c.add_capacitor("Co", "out", "gnd", p["Co"])
        c.add_resistor("Ro", "out", "gnd", p["Rload"])

    def set_frequency(self, freq):
        """Atualiza os gates para uma nova frequência"""
        dt = self.params["dead_time"]
        self.circuit.set_pwm("QH", freq, 0.5, 0.0, dt)
        self.circuit.set_pwm("QL", freq, 0.5, 0.5, dt)

    # =========================================================================
    # 1. SIMULAÇÃO TRANSIENTE (PARTIDA)
    # =========================================================================
    def run_transient(self, freq, cycles=150):
        print(f"\n--- 1. Rodando Transiente (Start-up) @ {freq/1e3:.1f} kHz ---")
        self.circuit.reset_state() # Garante estado zero
        self.set_frequency(freq)
        
        period = 1.0 / freq
        t_end = period * cycles
        
        # Roda simulação no tempo
        res = self.circuit.run(
            mode="transient", 
            dt=10e-9, 
            t_end=t_end,
            probes=["V(out)", "I(Lr)"]
        )
        
        # Plot
        t = np.array(res["t"]) * 1e3 # ms
        plt.figure(figsize=(10, 4))
        plt.subplot(1, 2, 1)
        plt.plot(t, res["signals"]["V(out)"], color='green')
        plt.title("Transiente de Tensão de Saída (Start-up)")
        plt.xlabel("Tempo (ms)"); plt.ylabel("Vout (V)")
        plt.grid(True)
        
        plt.subplot(1, 2, 2)
        # Zoom nos últimos ciclos
        zoom_idx = int(len(t)*0.9)
        plt.plot(t[zoom_idx:], res["signals"]["I(Lr)"][zoom_idx:], color='red')
        plt.title("Corrente Ressonante (Final)")
        plt.xlabel("Tempo (ms)"); plt.ylabel("ILr (A)")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    # =========================================================================
    # 2. SIMULAÇÃO SHOOTING (REGIME PERMANENTE DIRETO)
    # =========================================================================
    def run_shooting(self, freq):
        print(f"\n--- 2. Rodando Shooting Method @ {freq/1e3:.1f} kHz ---")
        self.set_frequency(freq)
        period = 1.0 / freq
        
        # A mágica acontece aqui: mode="shooting"
        t0 = time.time()
        res = self.circuit.run(
            mode="shooting",
            dt=10e-9,
            t_end=period*2, # Captura 2 ciclos para visualização
            period=period,
            probes=["V(out)", "I(Lr)", "V(tank)"]
        )
        print(f"Convergiu em {(time.time()-t0)*1000:.2f} ms!")
        
        # Plot
        t = np.array(res["t"]) * 1e6 # us
        plt.figure(figsize=(8, 4))
        plt.plot(t, res["signals"]["V(tank)"], label="V_Tank (Entrada)", alpha=0.6)
        plt.plot(t, np.array(res["signals"]["I(Lr)"])*10, label="I_Lr x10 (A)", color='red')
        plt.title(f"Regime Permanente (Shooting) @ {freq/1e3:.1f} kHz")
        plt.xlabel("Tempo (us)")
        plt.legend()
        plt.grid(True)
        plt.show()

    # =========================================================================
    # 3. SWEEP DE FREQUÊNCIA (CURVA DE GANHO)
    # =========================================================================
    def run_sweep(self):
        print("\n--- 3. Iniciando Frequency Sweep (50 pontos) ---")
        freqs = np.linspace(self.fr * 0.02, self.fr * 1.8, 100)
        gains = []
        
        self.circuit.set_params(max_iter=250, tol=1e-3, eps = 1) 

        for f in freqs:
            self.set_frequency(f)
            period = 1.0/f
            try:
                res = self.circuit.run(mode="shooting", dt= 10e-9, t_end=period, period=period, probes=["V(out)"])
                v_out = res["stats"]["V(out)"]["mean"]
                gains.append(v_out / self.params["Vin"])
            except:
                gains.append(0)
        
        plt.figure(figsize=(6, 4))
        plt.plot(freqs/1e3, gains, 'o-')
        plt.axvline(self.fr/1e3, color='r', linestyle='--', label='Fr')
        plt.title("Curva de Ganho LLC")
        plt.xlabel("Freq (kHz)"); plt.ylabel("Ganho (Vout/Vin)")
        plt.grid(True); plt.legend()
        plt.show()

    # =========================================================================
    # 4. ANÁLISE DE SENSIBILIDADE DA JACOBIANA (ESTABILIDADE)
    # =========================================================================
    def analyze_jacobian(self):
        print("\n--- 4. Análise de Sensibilidade da Jacobiana ---")
  
        test_freqs = [self.fr * 0.6, self.fr, self.fr * 1.4]
        
        print(f"{'Freq (kHz)':<15} | {'Cond. Number':<15} | {'Max Eigen (Abs)':<15} | {'Estabilidade'}")
        print("-" * 70)

        for f in test_freqs:
            self.set_frequency(f)
            period = 1.0/f
            
            res = self.circuit.run(
                mode="shooting", 
                dt = 10e-9, 
                period=period, 
                t_end=0, 
                export_matrices=True 
            )
            
            if "jacobian" in res:
                J = np.array(res["jacobian"]) 
                M = J + np.eye(J.shape[0])
                

                cond_num = np.linalg.cond(J)
                
                eig_vals = np.linalg.eigvals(M)
                max_eig = np.max(np.abs(eig_vals))
                
                status = "ESTÁVEL" if max_eig < 1.0 else "INSTÁVEL"
                
                print(f"{f/1e3:<15.1f} | {cond_num:<15.2e} | {max_eig:<15.4f} | {status}")
                
                # Plot dos autovalores no plano complexo (apenas para Fr)
                if f == self.fr:
                    plt.figure(figsize=(5, 5))
                    theta = np.linspace(0, 2*np.pi, 100)
                    plt.plot(np.cos(theta), np.sin(theta), 'k--', label='Círculo Unitário')
                    plt.scatter(eig_vals.real, eig_vals.imag, color='red', marker='x', s=100, label='Autovalores')
                    plt.title(f"Lugar das Raízes (Monodromia) @ {f/1e3:.1f} kHz")
                    plt.xlabel("Real"); plt.ylabel("Imag")
                    plt.grid(True); plt.legend()
                    plt.axis('equal')
                    plt.show()

# --- Execução ---
if __name__ == "__main__":
    sim = LLCFullAnalysis()
    
    # 1. Transiente
    sim.run_transient(freq=100e3)
    
    # 2. Shooting
    sim.run_shooting(freq=100e3)
    
    # 3. Sweep
    sim.run_sweep()
    
    # 4. Jacobiana
    sim.analyze_jacobian()