# Historia evolutiva de la ecolocalización y de la asociación con cuevas en murciélagos

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)

Este repositorio contiene los datos, el código fuente en R y los resultados íntegros generados para el análisis del Trabajo de Fin de Máster centrado en la historia evolutiva temprana de los murciélagos (Chiroptera). 

El objetivo principal de este trabajo es reconstruir la evolución de la ecolocación y la asociación con depósitos cavernícolas y kársticos, evaluando la robustez de las inferencias bajo dos marcos filogenéticos alternativos: la red morfológica masiva de Jones et al. (2024) y la topología molecular de Hand et al. (2023).

##  Estructura del repositorio

El análisis está diseñado para leer los datos de entrada y generar automáticamente la siguiente estructura de directorios de salida (`ASR_stochastic_mapping_outputs_publication/`):

*   **`data/`** (Archivos de entrada necesarios):
    *   `Data.xlsx`: Matriz de caracteres (Ecolocación, depósito en cueva, depósito kárstico) codificados bajo un modelo de incertidumbre (0, 0.5, 1) y rangos estratigráficos (FAD).
    *   `Hand_et_al_2023.tre`: Muestra posterior de árboles moleculares calibrados temporalmente.
    *   `Jones_et_al_2024.nex`: Árboles topológicos morfológicos sin datar.
*   **`scripts/`**:
    *   `A_Supplementary_ASR_stochastic_mapping_publication_ready.R`: Script principal ejecutable.
*   **`outputs/`** (Generados por el script):
    *   **`/Hand`** y **`/Jones`**: Resultados de la reconstrucción (AIC weights, `.csv` y `.rds` integrados).
    *   **`/Reference_ASR_trees`**: Reconstrucciones de estados ancestrales (ASR) proyectadas sobre el árbol MCC de Hand y el árbol de edades medias de Jones, incluyendo gráficos *pie chart* duales en formato PDF y PNG.
    *   **`/figures/08_probability_densitrees`**: Gráficos *densitree* de probabilidad enfrentados (Ostracoderm-style) generados con `ggtree` y `patchwork`.

##  Metodología Analítica

El código automatiza un flujo de trabajo filogenético comparativo completo:
1.  **Calibración temporal (Jones et al., 2024):** Generación estocástica de árboles datados en puntas utilizando la edad media y muestreos aleatorios dentro del rango estratigráfico (FAD_low - FAD_up) mediante el paquete `paleotree`.
2.  **Selección de modelos:** Comparación del ajuste empírico de modelos de transición evolutiva binaria: **ER** (Equal Rates), **SYM** (Symmetrical) y **ARD** (All-Rates-Different).
3.  **Mapeo Estocástico Integrado (SIMMAP):** Reconstrucción de estados ancestrales ponderada por el peso de Akaike (AICw) de los modelos para capturar la incertidumbre de las tasas de transición y la varianza topológica de los fósiles basales (*stem bats* como *Onychonycteris*, *Icaronycteris* y *Vielasia*).

##  Requisitos y dependencias

El script requiere **R (versión 4.0 o superior)**. El propio código incluye una función de autoinstalación para las dependencias faltantes. Los paquetes principales utilizados son:

**Análisis filogenético y comparativo:**
*   `ape`
*   `phytools`
*   `phangorn`
*   `paleotree`

**Manipulación de datos:**
*   `readxl`
*   `dplyr`
*   `tidyr`

**Visualización y gráficos:**
*   `ggplot2`
*   `ggtree` y `treeio` (vía Bioconductor)
*   `Ternary` (Gráficos de pesos AIC)
*   `patchwork` (Composición de figuras complejas)
*   `deeptime` (Escalas cronoestratigráficas)
*   `ggnewscale`, `scales`, `grid`


## Autoría

* **Dr. Humberto Gracián Ferrón**
  Universitat de València
  Departamento de Paleontología

* **Adriana Martínez Arroyo**
  Universitat de València
  Máster en Biodiversidad, Conservación y Evolución
