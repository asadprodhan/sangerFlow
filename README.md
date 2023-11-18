# **sangerFlow, an automated, reproducible, and scalable bioinformatics pipeline for analysing Sanger amplicon sequencing data leveraging Nextflow and Singularity** 

<br />

<br />

### **Asad Prodhan PhD** 


**https://asadprodhan.github.io/**


<br />


<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/sangerFlow/blob/main/sangerFlow.PNG"
  width=50% height=50%>
</p>
<p align = "center">
Figure 1: sangerFlow pipeline.
</p>



## **About the sangerFlow**


sangerFlow automatically analyses the forward and reverse reads from the Sanger amplicon sequencing data. The pipeline takes the fasta files as input and returns blastn hits i.e., species identifications for each amplicon. Therefore, the pipeline is automated and scalable.  Furthermore, the pipeline is written using the modern workflow manager, Nextflow; and Singularity containers. Therefore, it does not require software installation except Nextflow and Singularity, software subscription, or programming expertise from the end users.   
All these features make the pipeline ideal for large-scale Sanger amplicon sequencing data analysis and user-friendly. 


## **How to use the sangerFlow**


Follow the following steps to use sangerFlow.

### **Step 1: Install the required softwares**


-	Install conda

  
-	Create a conda environment and name it sangerFlow

  

```
conda create -n sangerFlow
```


-	Activate the conda environment sangerFlow


```
conda activate sangerFlow
```


-	Install Nextflow



```
conda install -c bioconda nextflow
```


-	Run the following command to make sure that Nextflow has been installed

  
```
nextflow -h
```


If you see the Nextflow options like Fig. 2, then the Nextflow has been installed



<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/sangerFlow/blob/main/Nextflow_options.PNG"
  >
</p>
<p align = "center">
Figure 2: Nextflow options.
</p>



-	Install Singularity


```
conda install -c conda-forge singularity
```


-	Run the following command to make sure that Singularity has been installed


```
singularity -h
```


If you see the Singularity options like Fig. 3, then the Singularity has been installed


<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/sangerFlow/blob/main/Singularity_options.PNG"
  >
</p>
<p align = "center">
Figure 3: Singularity options.
</p>


### **Step 2: Prepare a sample description file**


See Fig. 4. This is an example of a sample description file. It is a ‘tsv’ file format. 


-	First column is an Id for the sample

  
-	Second column is the forward (or read1) sequence file name

  
-	Third column is the reverse (or read2) sequence file name



<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/sangerFlow/blob/main/Sample_description_file.PNG"
  >
</p>
<p align = "center">
Figure 4: Sample description file.
</p>



### **Step 3: Download a blastn database from NCBI**


### **Step 4: Run sangerFlow**


-	Create a directory say ‘amplicon_analysis’

  
-	Transfer your Sanger amplicon sequencing data to ‘amplicon_analysis’ directory

  
-	Keep your sample description tsv file in the ‘amplicon_analysis’ directory

  
-	Run the following command to make sure the all the files are in UNIX format
  

```
dos2unix *
```


-	Run the following command to make sure that all the files are executable
  

```
chmod +x *
```


-	Run the following command to run sangerFlow

  

```
nextflow run asadprodhan/sangerFlow -r VERSION-NUMBER --db="/path/to/your/blastn_database"
```


> Collect the VERSION-NUMBER from the sangerFlow GitHub home page. It is located as shown in the red box in Fig. 5.



<br />
<p align="center">
  <img 
    src="https://github.com/asadprodhan/sangerFlow/blob/main/Version_number_location.PNG"
  >
</p>
<p align = "center">
Figure 5: sangerFlow version number location.
</p>


### **The End**
