library(Biostrings)
Archivo1 <- readDNAStringSet("~/Documents/Genomica/sequences/Arabidopsis_thaliana_nt.fasta")
Archivo2 <- readDNAStringSet("~/Documents/Genomica/sequences/Escherichia_coli_nt.fasta")
Archivo3 <- readDNAStringSet("~/Documents/Genomica/sequences/Saccharomyces_cerevisiae_nt.fasta")
Archivo4 <- readDNAStringSet("~/Documents/Genomica/sequences/Thermus_thermophilus_nt.fasta")
Archivo5 <- readDNAStringSet("~/Documents/Genomica/sequences/Homo_sapiens_nt.fasta")
chr1alf <- alphabetFrequency(Archivo1)
chr2alf <- alphabetFrequency(Archivo2)
chr3alf <- alphabetFrequency(Archivo3)
chr4alf <- alphabetFrequency(Archivo4)
chr5alf <- alphabetFrequency(Archivo5)
chr1gc <- sum(chr1alf[c("G", "C")]) / sum(chr1alf[c("A", "C", "G", "T")])
chr2gc <- sum(chr2alf[c("G", "C")]) / sum(chr2alf[c("A", "C", "G", "T")])
chr3gc <- sum(chr3alf[c("G", "C")]) / sum(chr3alf[c("A", "C", "G", "T")])
chr4gc <- sum(chr4alf[c("G", "C")]) / sum(chr4alf[c("A", "C", "G", "T")])
chr5gc <- sum(chr5alf[c("G", "C")]) / sum(chr5alf[c("A", "C", "G", "T")])


distriArchivo1 = (chr1alf[,2] + chr1alf[,3]) / (chr1alf[,1] + chr1alf[,2] + chr1alf[,3] + chr1alf[,4])
distriArchivo1
hist(distriArchivo1)


distriArchivo2 = (chr2alf[,2] + chr2alf[,3]) / (chr2alf[,1] + chr2alf[,2] + chr2alf[,3] + chr2alf[,4])
distriArchivo2
hist(distriArchivo2)


distriArchivo3 = (chr3alf[,2] + chr3alf[,3]) / (chr3alf[,1] + chr3alf[,2] + chr3alf[,3] + chr3alf[,4])
distriArchivo3
hist(distriArchivo3)



distriArchivo4 = (chr4alf[,2] + chr4alf[,3]) / (chr4alf[,1] + chr4alf[,2] + chr4alf[,3] + chr4alf[,4])
distriArchivo4
hist(distriArchivo4)

distriArchivo5 = (chr5alf[,2] + chr5alf[,3]) / (chr5alf[,1] + chr5alf[,2] + chr5alf[,3] + chr5alf[,4])
distriArchivo5
hist(distriArchivo5)

hist(distriArchivo1, distriArchivo2, distriArchivo3, distriArchivo4, distriArchivo5)

alf1 <- alphabetFrequency(chr1gc, collapse=TRUE)
alf2 <- alphabetFrequency(chr2gc, collapse=TRUE)
alf3 <- alphabetFrequency(chr3gc, collapse=TRUE)
alf4 <- alphabetFrequency(chr4gc, collapse=TRUE)
alf5 <- alphabetFrequency(chr5gc, collapse=TRUE)

df = data.frame(Organismo = c(rep("Archivo1", length(distriArchivo1)),
                rep("Archivo2", length(distriArchivo2))),
                GC = c(distriArchivo1, distriArchivo2))
library(ggplot2)
ggplot(df, aes(x = Organismo, y = GC)) + geom_boxplot()
