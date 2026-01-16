#pragma once
#include <cmath>
#include <limits>
#include <string>

class Waveform {
public:
    Waveform(const std::string& name, double t_start = 0.0)
        : name_(name), t_start_(t_start) {
    }

    void add(double t, double v) {
        if (t < t_start_) return;

        sum_ += v;
        sum2_ += v * v;
        min_ = std::min(min_, v);
        max_ = std::max(max_, v);
        count_++;
    }

    double mean() const {
        return count_ ? sum_ / count_ : 0.0;
    }

    double rms() const {
        return count_ ? std::sqrt(sum2_ / count_) : 0.0;
    }

    double min() const { return (count_ ? min_ : 0.0); }
    double max() const { return (count_ ? max_ : 0.0); }

    double ripple() const {
        return count_ ? (max_ - min_) : 0.0;
    }

    int samples() const { return count_; }
    const std::string& name() const { return name_; }

private:
    std::string name_;
    double t_start_;

    double sum_ = 0.0;
    double sum2_ = 0.0;
    double min_ = std::numeric_limits<double>::infinity();
    double max_ = -std::numeric_limits<double>::infinity();
    int count_ = 0;
};
