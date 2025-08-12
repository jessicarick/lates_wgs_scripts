# Scripts for Lates WGS and GBS data analysis

## Layout

The repository is split into seven main directories. Within each directory is a `README.md` file which summarizes the purpose of that directory; each section is also briefly described below. 

### **`experiments`** 
Where all of the experimental information lives, including any *executed* code. lives. This includes pipelines, scripts, and figure files. 
 * **`processing`**: Code used to *transform* the data into another type. 
 * **`analysis`**: Code to to *draw conclusions* from a data set.
 * **`exploratory`**: A sandbox that serves as a record of different approaches to transformation, interpretation, cleaning, or generation of data.
 * **`figures`**: Code used to generate figures.

### **`data`** 
Copies of the raw and processed data files that are small enough for hosting on GitHub. Larger data files can be found in the associated Zenodo repository (LINK). 

### **`miscellaneous`** 
Files that may not be code, but are important for reproducibility of your findings.
* **`protocols`**: Documentation of any relevant field, lab, etc. protocols. 

### **`software_module`** 
Custom code that is *not* executed directly, but is called from files in the `code` directory. 
