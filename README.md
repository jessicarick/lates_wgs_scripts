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


# License Information

<p xmlns:dct="http://purl.org/dc/terms/" xmlns:vcard="http://www.w3.org/2001/vcard-rdf/3.0#">
  <a rel="license"
     href="http://creativecommons.org/publicdomain/zero/1.0/">
    <img src="http://i.creativecommons.org/p/zero/1.0/88x31.png" style="border-style: none;" alt="CC0" />
  </a>
  <br />
  To the extent possible under law,
  <a rel="dct:publisher"
     href="github.com/gchure/reproducible_research">
    <span property="dct:title">Griffin Chure</span></a>
  has waived all copyright and related or neighboring rights to
  <span property="dct:title">A template for using git as a platform for reproducible scientific research</span>.
This work is published from:
<span property="vcard:Country" datatype="dct:ISO3166"
      content="US" about="github.com/gchure/reproducible_research">
  United States</span>.
</p>
