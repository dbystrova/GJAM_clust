
gjamSpec2Trait <- function(pbys, sbyt, tTypes){
  
  # plotBySpec  - n by S numeric matrix
  # specByTrait - S by M data.frame
  # traitTypes  - data types for traits
  # FC can be factors that will be categorical
  
  if(!identical(colnames(pbys),rownames(sbyt))){
    stop( 'colnames(pbys) must match rownames(sbyt)' )
  }
  if(length(tTypes) != ncol(sbyt)){
    stop( 'length(tTypes) must equal ncol(sbyt)' )
  }
  
  .spec2Trait(pbys, sbyt, tTypes)
}
                