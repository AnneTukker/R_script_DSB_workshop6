library(tidyverse) #load tidyverse package
#Reading in 3 data files and attaching them to objects so it is easier to refer to in the rest of the code:
genotypes <- read_tsv("~/OneDrive - University of East Anglia/MSc molecular medicine/Datascience/workshop6/GWAS_data/Chr12_Genotypes.tsv") 
map <- read_tsv("~/OneDrive - University of East Anglia/MSc molecular medicine/Datascience/workshop6/GWAS_data/Chr12_Map.tsv") 
phenotypes <- read_tsv("~/OneDrive - University of East Anglia/MSc molecular medicine/Datascience/workshop6/GWAS_data/Phenotype_Height.tsv") 
#CHECK column names in the files:
colnames(genotypes) 
colnames(map)
colnames(phenotypes)
#CHECK number of rows in the files:
nrow(genotypes) 
nrow(map)
nrow(phenotypes)
#CHECK number of columns in the files:
ncol(genotypes) 
ncol(map)
ncol(phenotypes)
#MAKING a histogram of variable height:
ggplot(phenotypes, aes(x = Height)) + geom_histogram() #geom_histogram() automatically chooses a bin width, you can specify this yourself by adding a binwidth argument
#MAKING a histogram to check if sites are evenly distributed along the chromosome:
ggplot(map, aes(x=Position)) + geom_histogram()
#DETERMINE allele frequency:
one_geno = as.numeric(genotypes[1,]) #leaving the column specification blank to give me the whole row, I make it a vector to calculate the allele frequency
one_geno 
mean(one_geno) #there are some N/A values -> can't average that so the output is N/A -> we need to exclude N/A values
mean(one_geno, na.rm = TRUE) #the added argument ensures it excludes N/A values
mean(one_geno, na.rm = TRUE)/2 #/2 because we code the genotype as '0' or '2'
#APPLYING function to every single row:
frequencies = apply(genotypes,1,mean,na.rm=TRUE)/2 #1 to specify we want to apply it over every signle row (2 would be every single columm)
head(frequencies) #to view the top of the frequencies value
plot_data = tibble(frequencies) #creating a table so we can make a plot
ggplot(plot_data, aes(x=frequencies)) + geom_histogram(binwidth = 0.05) #to create a histrogram with binwidth 0.05
ggplot(plot_data, aes(x=frequencies)) + geom_histogram(binwidth = 0.02) #data is skewed, I adjusted the binwidth to create more detail in the histogram
#the "plot_data =" step (line 31) can also be incorporated in the code in line 32 so you code it all at once
#PLOTTING association between genotype and phenotype:
phenotypes$one_geno = one_geno #adding a column (called one_geno) to the phenotypes table
ggplot(phenotypes, aes(x=one_geno,y=Height)) + geom_jitter(width=0.2) #geom_jitter to visualize all the datapoints in the plot
ggplot(phenotypes, aes(x=one_geno,y=Height)) + geom_jitter(width=0.2) + geom_smooth(method = 'lm') #geom_smooth(method = 'lm') to add a line in the plot
#STATISTICAL association:
corTestResult = cor.test(phenotypes$Height, phenotypes$one_geno)
corTestResult
corTestResult$p.value
#AUTOMATING statistical association:
association_test = function(genotype, phenotype){
  genotype = as.numeric(unlist(genotype))
  phenotype = as.numeric(unlist(phenotype))
  
  cor_test_result = cor.test(genotype, phenotype) 
  cor_test_result$p.value
}
apply(genotypes,1,association_test,phenotype=phenotypes$Height) #automate the repeat of the function (association_test) by row

map$pvalues = apply(genotypes,1,association_test,phenotype=phenotypes$Height) #add the results (pvalues) of the automation in a new column in the map table
ggplot(map, aes(x=Position, y=-log10(pvalues)))+geom_point() #create a plot with position on the X axis and the pvalues (-log10) on the Y axis
