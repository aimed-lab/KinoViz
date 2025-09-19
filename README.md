KinoViz: Interactive Kinome Profiling Analysis

    KinoViz is a user-friendly, interactive web application designed for the analysis and visualization of high-throughput kinome profiling data. Built with R and Shiny, it empowers cancer researchers and biologists to explore complex peptide phosphorylation patterns without needing advanced programming skills.

This tool bridges the gap between raw experimental data and actionable biological insights by providing a suite of powerful visualization and statistical modules.

![alt text](https://via.placeholder.com/800x450.png?text=KinoViz+Main+Dashboard+Screenshot)


(Placeholder: Replace with an actual screenshot of the application's welcome screen)
✨ Key Features

    Interactive Data Exploration: Upload, view, and clean your dataset in a responsive data grid.

    Dynamic Filtering: Filter data globally by Sample, Array, Cycle, and Exposure Time to focus on specific experimental conditions.

    Slope Correction: Apply a linear correction to signal intensities to account for variations across multiple exposure times.

    Top Peptide Identification: Automatically identify and visualize the top N peptides based on maximum signal intensity, or select specific peptides for analysis.

    Kinetic Analysis: Visualize peptide phosphorylation kinetics over time (Cycle) or exposure time.

    Comparative Analysis: Compare signal profiles of multiple selected peptides simultaneously.

    Clustered Heatmaps: Generate globally clustered heatmaps to visualize signal patterns across samples and peptides.

    Dimensionality Reduction: Perform PCA and UMAP analysis to visualize sample clustering, with options to group by Sample or Array.

    Network Analysis: Construct and visualize a peptide correlation network to identify clusters of co-regulated peptides. The network includes an interactive legend for edge weights.

    Statistical Testing: A dedicated module to perform one-way ANOVA and post-hoc Tukey's HSD tests to statistically compare peptide signals between different sample groups.

    Publication-Ready Plots: Export all visualizations as high-resolution PNG images suitable for publications, with user-configurable DPI settings.

🛠️ Technical Stack

    Backend: R

    Web Framework: Shiny

    Core Packages:

        plotly for interactive visualizations

        ggplot2 for static plots

        dplyr & tidyr for data manipulation

        igraph for network analysis

        umap for dimensionality reduction

        DT for interactive tables

🚀 Getting Started

Follow these instructions to set up and run KinoViz on your local machine.
Prerequisites

    R (version 4.0 or newer)

    RStudio IDE (recommended)

Installation

    Clone the Repository:
    code Bash

IGNORE_WHEN_COPYING_START
IGNORE_WHEN_COPYING_END

    
git clone https://github.com/your-username/KinoViz.git
cd KinoViz

  

Install R Packages:
Open the app.R file in RStudio. The script will automatically detect the required packages. Alternatively, you can run the following command in your R console to install all dependencies at once:
code R

    IGNORE_WHEN_COPYING_START
    IGNORE_WHEN_COPYING_END

        
    install.packages(c(
      "shiny", "shinydashboard", "shinythemes", "shinyWidgets", "shinyjs",
      "readr", "data.table", "dplyr", "tidyr", "stringr", "reshape2",
      "ggplot2", "scales", "plotly", "heatmaply", "ggridges", "RColorBrewer",
      "viridis", "igraph", "umap", "jsonlite", "aggrid", "DT"
    ))

      

    Prepare Data:
    The application is configured to load a default dataset named clean_kinome.csv. Ensure this file is present in the same directory as app.R. If you wish to use your own data, see the Data Format section below.

Running the Application

    Open the app.R file in RStudio.

    Click the "Run App" button at the top-right of the script editor.

Alternatively, you can run the app from the R console:
code R
IGNORE_WHEN_COPYING_START
IGNORE_WHEN_COPYING_END

    
# Make sure your working directory is set to the project folder
shiny::runApp()

  

📋 Usage Workflow

    Load Data: Use the sidebar to upload your own CSV file or click "Load & Process Data" to use the default clean_kinome.csv.

    Apply Filters: Once the data is loaded, filtering options will appear in the sidebar. Select the Sample(s), Array(s), Cycle, and Exposure Time you wish to analyze. These filters apply globally to most tabs.

    Explore Analysis Tabs: Navigate through the different tabs to perform analyses:

        Explore & Clean: View the raw data and apply slope correction.

        Top Peptides: Find high-signal peptides.

        Dim Reduction: Run PCA/UMAP to see how your samples cluster.

        Statistical Analysis: Perform ANOVA to compare your sample groups.

        ...and more.

    Export Results:

        For interactive plots, use the camera icon in the plot's mode bar to save a high-resolution PNG. The resolution can be adjusted using the "PNG Scale Factor" in the sidebar.

        For static plots (like the one in "Top Peptides"), use the dedicated download button to save a high-resolution image.

        Tables and processed data can also be downloaded using their respective download buttons.

📄 Data Format

To use your own data with KinoViz, your CSV file must adhere to the following format:

    The file must have a header row.

    The first four columns must be, in order:

        Sample: The name of the treatment, cell line, or sample group (e.g., "Control", "DrugA").

        Array: The identifier for the technical replicate or array (e.g., "A1", "A2").

        Cycle: The measurement cycle number (numeric).

        Exposure_time: The camera exposure time (numeric).

    All subsequent columns should contain the signal intensity data for each peptide, with the peptide ID as the column header.

Sample	Array	Cycle	Exposure_time	Peptide_A	Peptide_B	...
Control	A1	1	100	15023.4	8765.1	...
Control	A1	2	100	18045.2	9123.5	...
DrugA	A1	1	100	45011.9	10234.0	...
...	...	...	...	...	...	...
📜 Citation

If you use KinoViz in your research, please cite:

    Saghapour, E., Anderson, J. C., Chen, J., & Willey, C. D. (2025). KinoViz: A User-Friendly Web Application for High-Throughput Kinome Profiling Analysis and Visualization in Cancer Research. [Provide URL or Journal Reference when available].

🤝 Contributing

Contributions are welcome! If you have suggestions for improvements or encounter a bug, please open an issue.

To contribute code:

    Fork the repository.

    Create a new branch (git checkout -b feature/your-feature).

    Commit your changes (git commit -m 'Add some feature').

    Push to the branch (git push origin feature/your-feature).

    Open a Pull Request.

📄 License

This project is licensed under the MIT License. See the LICENSE file for details.
📧 Contact

    Principal Investigator: Christopher D. Willey, MD, PhD - cwilley@uabmc.edu

    Lead Developer: Ehsan Saghapour
