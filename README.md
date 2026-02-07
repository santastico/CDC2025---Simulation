# Controllers Synthesis for Max-Plus Linear Systems under State Constraints

[![CDC 2025](https://img.shields.io/badge/CDC-2025-blue.svg)](https://cdc2025.ieeecss.org/)
[![Max-Plus Algebra](https://img.shields.io/badge/Max--Plus-Algebra-green.svg)]()
[![License](https://img.shields.io/badge/license-MIT-orange.svg)]()

This repository contains the presentation and simulation tools for the paper **"Controllers Synthesis for Max-Plus Linear Systems under State Constraints"** presented at the 64th IEEE Conference on Decision and Control (CDC 2025).

### Authors
- **G. S. Pereira** - LARIS, Université d'Angers, France
- **L. Hardouin** - LARIS, Université d'Angers, France
- **B. Cottenceau** - LARIS, Université d'Angers, France
- **C. A. Maia** - LOPAC, DEE-EE, UFMG, Brazil

---

## 📋 Overview

This work addresses the synthesis of neutral controllers for **Max-Plus Linear Systems** (modeling Timed Event Graphs - TEGs) under state constraints. The proposed methodology ensures that:

- The controlled system preserves optimal input/output performance
- State trajectories satisfy two-sided equations $Nx = Mx$
- Controllers are causal and realizable
- The control law $u = P(v \oplus Fy)$ combines prefilter ($P$) and feedback ($F$) components

The max-plus algebraic framework provides powerful tools for modeling and analyzing discrete event systems with synchronization and delay constraints.

---

## 📁 Repository Structure

```text
CDC2025---Simulation/
│
├── Apresentacao_CDC2025.pdf          # Conference presentation slides
├── Obtencao_FandP.cpp                # C++ program to compute controllers
├── Simulacao_HTML_CDC.html           # Main HTML visualization interface
├── Simulacao_SVG.svg                 # SVG graphics for TEG visualization
├── jquery-3.7.1.slim.min.js          # jQuery library for interactivity
└── README.md                         # This file
