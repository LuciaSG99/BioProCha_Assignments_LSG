class InteractionNetwork
#contains the members of each network. 

  attr_accessor :gene_original #The object of the first gene used to search for interactions
  attr_accessor :network 
  attr_accessor :all_genes #hash_interactions
  attr_accessor :go_annotations
  attr_accessor :kegg_annotations
  require 'rest-client'
  require 'json'
  
  def initialize(params)
    @gene_original = params.fetch(:gene_original, nil)
    @all_genes = params.fetch(:all_genes, Hash.new)
    @network = Array.new
    recursive_function(@gene_original)
    @go_annotations = Hash.new
    @kegg_annotations = Hash.new
    obtain_kegg_network
    obtain_go_network
    
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
  
  def fetch(url, headers = {accept: "*/*"}, user = "", pass="")
      response = RestClient::Request.execute({
        method: :get,
        url: url.to_s,
        user: user,
        password: pass,
        headers: headers})
      return response
      
      rescue RestClient::ExceptionWithResponse => e
        $stderr.puts e.inspect
        response = false
        return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
      rescue RestClient::Exception => e
        $stderr.puts e.inspect
        response = false
        return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
      rescue Exception => e
        $stderr.puts e.inspect
        response = false
        return response  # now we are returning 'False', and we will check that with an \"if\" statement in our main code
  end
  
  
  
  def obtain_kegg_network
    @network.each{|gene_id| result_kegg_patways = fetch("http://togows.org/entry/kegg-genes/ath:#{gene_id.agi_locus}/pathways.json")
        kegg_patways_json = JSON.parse(result_go_patways) #transform the results into JSON format
        kegg_pathways_gene = kegg_patways_json[0]  #store the hash with the KEGG ID and the patways' name
        unless @kegg_annotations.has_key?
            #code
        end
        
        kegg_pathways_gene.each {|key,value| unless @kegg_annotations.has_key?(key)
            @kegg_annotations[key] = value
        end}}
  end
  
  def obtain_go_network
    @network.each{|gene_id| result_go_patways = fetch("http://togows.org/entry/ebi-uniprot/#{gene_id.agi_locus}/dr.json")
    go_pathways_json = JSON.parse(result_go_patways) #transform the results into JSON format
    #select only the biological process:
    go_biological_process = go_pathways_json[0]["GO"].select!{|array| array.any?{|element| element =~ /P:/}}
    go_biological_process.each{|result| @go_annotations[result[0]] = result[1]}}
  end
  
  
#  def obtain_gene_kegg(gene_id)
#    result_kegg_patways = fetch("http://togows.org/entry/kegg-genes/ath:#{gene_id}/pathways.json")
#    kegg_patways_json = JSON.parse(result_go_patways) #transform the results into JSON format
#    kegg_pathways_gene = kegg_patways_json[0] #{"ath00270"=>"Cysteine", "ath01100"=>"Metabolic pathways"}
#    
#    return 
#  end
#  
#  def obtain_gene_go(gene_id)
#    result_go_patways = fetch("http://togows.org/entry/ebi-uniprot/#{gene_id}/dr.json")
#    go_patways_json = JSON.parse(result_go_patways) #transform the results into JSON format
#    return go_pathways_json
#  end
#  
end
