#' getConversion Method
#'
#' @param object CoSIAn. 
#'
#' @return 
#' @export
#'
#' @examples
setMethod("getConversion", signature(object = "CoSIAn"), function(object) { # user's input of the function
  #Set each part of the object that this method uses into their own variable that will be used inside the code
    input_species<-object@i_species
    input_id<-object@i_id
    input<-object@gene_set
    output_ids<-object@o_ids
    output_species<-object@o_species
    tool<-object@mapping_tool
    ortholog_database<-object@ortholog_database
    Filter_I_Species <- switch(input_species,
                               homo_sapiens = {
                                 species_data<-homo_sapiens(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               mus_musculus = {
                                 species_data<-mus_musculus(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               r_norvegicus = {
                                 species_data<-r_norvegicus(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               danio_rerio = {
                                 species_data<-danio_rerio(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               c_elegans = {
                                 species_data<-c_elegans(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               d_melanogaster = {
                                 species_data<-d_melanogaster(input_id,input,output_ids,output_species, tool, ortholog_database)
                               },
                               stop("Error: Invalid i_species in CoSIAn Object. Make sure the species in the i_species slot is an avalible model organism and is in the correct format.")
                               
    )
    return(species_data)
})




