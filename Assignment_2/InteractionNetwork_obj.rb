class InteractionNetwork
#contains the members of each network. 

  attr_accessor :gene_original #The object of the first gene used to search for interactions
  attr_accessor :network 
  attr_accessor :all_genes #hash_interactions
  
  def initialize(params)
    @gene_original = params.fetch(:gene_original, nil)
    @all_genes = params.fetch(:all_genes, Hash.new)
    @network = Array.new
    recursive_function(@gene_original)
  end
      
  
  def recursive_function(gene_obj)
    #puts @network.length
    if @network.include? (gene_obj)
        return    
    end         
    @network << gene_obj #mete el objeto gen en el atributo network que es una lista
    gene_obj.interactions.each {|gene_int| #para cada interaccion (que es una string)
        if @all_genes.has_key? (gene_int.downcase) #si ese gen es una key del diccionario con los genes objetos (es decir, si es un objeto Gen)
            int_obj = @all_genes[gene_int.downcase] #dame ese objeto con sus respectivos atributos Gen1.2.(objeto)           
            recursive_function(int_obj)
        
        end}      
      
  end
end
