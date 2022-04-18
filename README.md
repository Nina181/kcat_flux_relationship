# kcat_flux_relationship

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#Packages">Packages</a></li>
        <li><a href="#First use">First use</a></li>
      </ul>
    </li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

Source code of the bachelor thesis "Analyzing the Relationship Between Enzyme Turnover Numbers and Metabolic Fluxes Predicted by Flux Balance Analyses".

<!-- GETTING STARTED -->
## Getting Started

To set up this project locally follow these steps. 
Python version 3.8 was used for this project.

### Packages

Used packages:
* pandas==1.14
* cobra==0.22.1
* bio==1.3.4
* ete3==3.1.2
* matplotlib==3.1.1
* scipy==1.4.1
* requests==2.25.1
  

### First use

1. Clone the repo
   ```sh
   git clone https://github.com/Nina181/kcat_flux_relationship
   ```
2. Install packages
   ```sh
   pip install <package_name>
   ```
3. Update NCBI Database (only necessary for first use)
   ```js
   uncomment line 93 in method add_lineage() of the class kcat_flux_mapping.mapping
   ```
4. Enter your email address necessary for accessing the NCBI Database
   ```js
   in method add_lineage() line 94 of the class kcat_flux_mapping.mapping
   ```

<p align="right">(<a href="#top">back to top</a>)</p>
