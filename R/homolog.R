homolog<-function(entrez_data, species_number, ortholog_database){
  if (ortholog_database=="HomoloGene"){
    myGenes<- as.character(entrez_data)
    species_one<-entrez_data
    data<- annotationTools::getHOMOLOG(myGenes, species_number,homologene, noIDsymbol = NA)
    data<-as.list(data)
    species_two<- as.character(data)
    homologs<- as.data.frame(species_two)
    homologs<-merge(data.frame(homologs, row.names=NULL), data.frame(species_one, row.names=NULL),
                    by = 0, all = TRUE)[-1]
    manipulated_homologs<- data.frame(homologs)
    x<-1:nrow(manipulated_homologs)
    for (i in seq_along(x)) {
      ##print(i)
      char<-manipulated_homologs[i,1]
      ##print(char)
      if ((char=="NA")==TRUE){
        manipulated_homologs[i,1]<-NA
      }
    }
    species_two_entrez_ID<-na.omit(manipulated_homologs)
    return(species_two_entrez_ID)
  }
  if (ortholog_database=="NCBIOrthoAnnotationPipe"){
    myGenes<- as.vector(entrez_data)
    species_one<-entrez_data
    data<- annotationTools::getHOMOLOG(myGenes, species_number,NCBIOrtho, tableType="gene_orthologs")
    data<-as.list(data)
    species_two<- as.character(data)
    homologs<- as.data.frame(species_two)
    homologs<-merge(data.frame(homologs, row.names=NULL), data.frame(species_one, row.names=NULL),
                    by = 0, all = TRUE)[-1]
    manipulated_homologs<- data.frame(homologs)
    x<-1:nrow(manipulated_homologs)
    for (i in seq_along(x)) {
      ##print(i)
      char<-manipulated_homologs[i,1]
      ##print(char)
      if ((char=="NA")==TRUE){
        manipulated_homologs[i,1]<-NA
      }
    }
    species_two_entrez_ID<-na.omit(manipulated_homologs)
    return(species_two_entrez_ID)
  }