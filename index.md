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

Insert full path to csv file with the input data, choose MDL or NML using the toggle button, then specify the number of random initializations and the maximum number of clusters.

![](https://i.imgur.com/qjrg9Iz.png)

![](https://i.imgur.com/AWKv6ur.png)

![](https://i.imgur.com/vBIcoyZ.png)
