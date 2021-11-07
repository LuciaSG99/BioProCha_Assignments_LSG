class InteractionNetwork
#contains the members of each network. 

  attr_accessor :gene_original #The object of the first gene used to search for interactions
  attr_accessor :network 
  attr_accessor :all_genes #hash_interactions
  
  def initialize(params)
    @gene_original = params.fetch(:gene_original, nil)
    @all_genes = params.fetch(:all_genes, Hash.new)
    @network = recursive_function(@gene_original,Array.new,@all_genes,contador = 0)
  end
      
  
  def recursive_function(gene_obj,all_genes)  #contador
    if @network.include? (gene_obj.agi_locus)
        return    
    else
      #if contador > 3
      #    return 
      #else
      if gene.type == Gene        
        @network << gene_obj #.agi_locus
        interactions = gene_obj.interactions #es una lista con las interacciones del gene
        #puts interactions
      else 
        if interactions.empty?
          return 
        else  
          interactions.each {|gene_int| 
            if all_genes.has_key? (gene_int.downcase) 
              int_obj = all_genes[gene_int.downcase]
              #puts int_obj
              puts 
              @network << int_obj #.agi_locus
              
              if int_obj.interactions.empty?
                next
              else
                  recursive_function(int_obj,all_genes) #,contador+1)
              end
            else          
              @network << gene_int.downcase
              next           
              
            end}
          end
      end
    end 
  end
end 