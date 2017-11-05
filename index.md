---
layout: default
---

#### Program description

pkclusteringmdl is a Java implementation for clustering of pharmacokinetics responses with model selection based on the Minimum Description Length principle or the Normalized Maximum Likelihood codelength. It takes a CSV file and performs clustering using one of the model selection techniques according to a user-defined number of random initializations and maximum number of output clusters. The result is shown in a new window on-screen but can be saved in a text file Result.txt.

#### Usage

To use the program, execute the jar file:

```
$ java -jar pkclusteringmdl.jar
```

Upon execution of the program, an interface window as seen on the image below will appear. In this window, the user may insert the path to a CSV file with the pharmacokinetic data to be used. There is a toggle button for if the user decides to use NML. If this button is not pressed, MDL will be used by default.

![](https://i.imgur.com/qjrg9Iz.png)

If the file was read correctly, a new window will appear as seen below. In this window, it is possible to set the maximum number of clusters that can be given as the output of the program, as well as the number of random initializations of the EM algorithm to be used for each of the possible numbers of clusters. To execute the algorithm, press Run.

![](https://i.imgur.com/AWKv6ur.png)

The program will then run the algorithm for the data and parameters provided using a parallel implementation. When the execution is finished, a new window with the results will appear as seen below. Pressing the Save button will save the results in a text file named "Result.txt" in the same directory as the jar file. Rename this file in order to save multiple result files.

![](https://i.imgur.com/vBIcoyZ.png)
