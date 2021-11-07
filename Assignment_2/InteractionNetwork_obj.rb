class InteractionNetwork
#contains the members of each network. 

  attr_accessor :gene_original #The object of the first gene used to search for interactions
  attr_accessor :network ["auadtud"]
  attr_accessor :all_genes #hash_interactions
  #attr_accesor :level_1_genes
  #attr_accesor :level_2_genes
  #attr_accesor :level_3_genes
  #
  def initialize(params)
    @gene_original = params.fetch(:gene_original, nil)
    @all_genes = params.fetch(:all_genes, Hash.new)
    @network = recursive_function(@gene_original,Array.new,@all_genes,contador = 0)
  end
      
  
  def recursive_function(gene_obj,network,all_genes,contador)  
    if network.include? (gene_obj.agi_locus)
        return network   
    else
      if contador > 3
          return network
      else
        network << gene_obj.agi_locus
        interactions = gene_obj.interactions #es una lista con las interacciones del gene
        #puts interactions
        if interactions.empty?
          return network
        else  
          interactions.each {|gene_int| #gene_int = str 
            #puts                                          Gen1 --> int : Gen1.2 --> Gen1.2.3 --> Gen1.2.3.4
            #puts gene_int
            if all_genes.has_key? (gene_int.downcase) 
              int_obj = all_genes[gene_int.downcase]
              #puts int_obj
              puts 
              network << int_obj.agi_locus
              
              if int_obj.interactions.empty?
                next
              else
                  recursive_function(int_obj,network, all_genes,contador+1)
              end
            else          
              network << gene_int.downcase
              next           
              
            end}
          return network
        end
      end
    end 
  end
end 