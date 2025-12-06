# ES313 — Simulation & Modeling in Julia

This repository contains the modules, scripts, and tools developed for the **ES313 Simulation Project**, focusing on orbital mechanics, satellite constellation modeling, visibility analysis, and Earth-coverage computation using Julia.  
The objective is to provide a clean, modular, and reusable framework for simulating **LEO constellations**, evaluating coverage performance, and generating 2D/3D visualisations.

---

## 📁 Repository Structure

### **📦 Core Julia Module**
- **`SatsLEO.jl`**  
  Implements:
  - ECI/ECEF coordinate transforms  
  - Rotation matrices  
  - Walker-Delta constellation generation  
  - Custom constellation layouts (via `myconstellation`)  
  - Line-of-sight visibility checks  
  - Instantaneous and average coverage computation  
  - 2D/3D visualisation tools  

### **📁 Notebooks / Scripts**
- Constellation generation & 3D visualisation  
- Instantaneous and mean coverage calculations  
- Global coverage heatmaps (latitude–longitude grids)  

---

## 🌐 Code Documentation Language

**All code comments and docstrings inside the Julia source files are written in French.**  
The README remains in **English** for broader accessibility.

---

## 📦 Dependencies

Minimal Julia dependencies:

- `LinearAlgebra`  
- `Plots.jl`

Multithreading (optional): `Threads.@threads`
Logging (optinal): `Logging`

---

## 🛰️ Example Usage (LEO Constellations)

```julia
using .SatsLEO

P, S, F = 4, 12, 1        # Walker-Delta parameters
i_deg = 50                # Orbital inclination (degrees)
h_km = 800                # Satellite altitude above Earth (km)
a = Re + h_km * 1e3       # Semi-major axis

sats = walker_delta(P, S, F, i_deg, a)     # Generate the constellation

mean_coverage_fraction(sats, -i_deg, i_deg, 10)
# → Computes average coverage between ±inclination with a 10° elevation mask
```

--- 

## 📜 License

**MIT License**

---

## 🤝 Contributors

Developed by **Aurelien Habets**, with assistance from **OpenAI ChatGPT 5.**