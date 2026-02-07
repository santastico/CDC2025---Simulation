## 📄 File Descriptions

### 1. **Apresentacao_CDC2025.pdf**
**Conference Presentation Slides**

The complete presentation delivered at CDC 2025, containing:
- Introduction to Timed Event Graphs (TEGs) and max-plus algebra
- Problem formulation for controller synthesis under state constraints
- Theoretical framework: residuation theory, Kleene star operator, γ-transform
- Methodology for computing optimal prefilter $P_{opt}^+$ and feedback $F_{opt}^+$
- Numerical example with TEG model and constraint specifications
- Visualization of the controlled system

**Key Concepts:**
- Max-plus semiring: $(\mathbb{Z}_{max}, \oplus, \otimes)$ where $a \oplus b = \max(a,b)$ and $a \otimes b = a + b$
- Dynamic model: $x(k) = Ax(k-1) \oplus Bu(k)$, $y(k) = Cx(k)$
- Neutral control: Transfer matrix preserved $H = (CA^*BPF)^*CA^*BP \preceq CA^*B$

### 2. **Obtencao_FandP.cpp**
**Controller Computation Program**

C++ implementation for computing the optimal prefilter and feedback controllers.

**Functionality:**
- Computes initial prefilter parameter: $P_{opt} = (CA^*B) \backslash^{\circ} (CA^*B)$
- Applies iterative algorithm to find $P_{opt}^+$ using the mapping:
  $$f(P) = Pr_+(P) \wedge (NA^*B) \backslash^{\circ} (MA^*BP) \wedge (MA^*B) \backslash^{\circ} (NA^*BP)$$
- Calculates optimal feedback: $F_{opt}^+ = Pr_+(P_{opt} \backslash^{\circ} P_{opt} /^{\circ} (CA^*BP_{opt}))$
- Handles matrix operations in max-plus algebra (residuation, Kleene star, causality projection)

**Input:** System matrices A, B, C and constraint matrices N, M  
**Output:** Controller matrices $P_{opt}^+$ and $F_{opt}^+$ in series form

**Example Output (from presentation):**
```text
P_opt^+ = [ (12γ⁴ ⊕ 20γ⁵)r*    (8γ⁰ ⊕ 14γ²)r* ]
          [ (4γ⁴ ⊕ 12γ⁵)r*     (0γ⁰ ⊕ 6γ²)r*  ]

F_opt^+ = [ (3γ¹ ⊕ 9γ³)r* ]
          [ (1γ³ ⊕ 7γ⁵)r* ]

where r = (14γ³)
```

### 3. **Simulacao_HTML_CDC.html**
**Interactive HTML Simulation Interface**

Web-based visualization tool for simulating the controlled TEG system.

**Features:**
- Interactive TEG graph visualization with transitions and places
- Real-time simulation of firing sequences
- Visual representation of token flow and timing constraints
- Display of state trajectories $x(k)$ under control
- Verification of constraints $Nx = Mx$
- Comparison between controlled and uncontrolled behavior

**Usage:**
1. Open `Simulacao_HTML_CDC.html` in a web browser
2. Input reference signal $v(k)$
3. Observe controlled system response $y(k)$
4. Verify constraint satisfaction visually

**Technologies:** HTML5, CSS3, JavaScript, SVG manipulation

### 4. **Simulacao_SVG.svg**
**TEG Graphical Representation**

Scalable Vector Graphics file containing the visual model of the Timed Event Graph.

**Content:**
- Transitions (circles) with timing information
- Places (bars) with token markings
- Weighted arcs representing synchronization and delays
- Controller blocks (P and F) integrated into the TEG structure
- Reference input $v$ and output $y$ connections

**Integration:** Loaded by `Simulacao_HTML_CDC.html` for dynamic animation

**Example System (from presentation):**
- 4 internal transitions: $x_1, x_2, x_3, x_4$
- 2 input transitions: $u_1, u_2$
- 1 output transition: $y_1$
- Constraints: $x_2 - x_1 \leq 6$ and $x_2 - x_3 \leq 12$

### 5. **jquery-3.7.1.slim.min.js**
**jQuery Library**

Minified jQuery library (version 3.7.1 slim) used for DOM manipulation and event handling in the HTML simulation.

**Purpose:**
- Simplifies JavaScript interactions with HTML elements
- Handles user input events (button clicks, parameter changes)
- Manages SVG animation and state updates
- Provides cross-browser compatibility

**Note:** Slim version excludes AJAX and effects modules for reduced file size.

---

## 🚀 Getting Started

### Prerequisites
- **For Controller Computation:** C++ compiler (g++, clang++)
- **For Visualization:** Modern web browser (Chrome, Firefox, Edge, Safari)

### Running the Controller Computation

```bash
# Compile the C++ program
g++ -std=c++11 -o controller Obtencao_FandP.cpp

# Run the program
./controller
```

The program will output the optimal controller matrices $P_{opt}^+$ and $F_{opt}^+$ in max-plus series notation.

### Running the Simulation

1. Ensure all files are in the same directory
2. Open `Simulacao_HTML_CDC.html` in a web browser
3. The TEG visualization will load from `Simulacao_SVG.svg`
4. Interact with the simulation controls to observe system behavior

**No server required** - runs entirely client-side!

---

## 🤝 Acknowledgments

- **Departamento de Engenharia Elétrica, Escola de Engenharia, Universidade Federal de Minas Gerais (UFMG)**
- **Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq)**
- **École polytechnique de l'université d'Angers (Polytech Angers)**

---

## 📧 Contact

For questions, collaboration opportunities, or issues with the simulation, please contact:
- **G. S. Pereira** - [santospereiragsp@gmail.com](mailto:santospereiragsp@gmail.com)

---

## 📝 License

This project is available for free under an open-source license. Please cite our CDC 2025 paper if you use this work in academic research.

---

## 🔗 Links

- **GitHub Repository:** [https://github.com/santastico/CDC2025---Simulation](https://github.com/santastico/CDC2025---Simulation)
- **CDC 2025 Conference:** [https://cdc2025.ieeecss.org/](https://cdc2025.ieeecss.org/)

---

**Last updated:** February 2026# Controllers Synthesis for Max-Plus Linear Systems under State Constraints

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
