## main code
file = ARGV[2]
require 'rest-client'
require '/home/osboxes/BioinformaticsCourseGitHub/BioProCha_Assignments_LSG/Assignment_2/InteractionNetwork_obj.rb'
require '/home/osboxes/BioinformaticsCourseGitHub/BioProCha_Assignments_LSG/Assignment_2/gene_obj.rb'



#saving the genes of the list in an array:
def read_file_genes(file)
  genes_list = Array.new
  lines = IO.readlines(file)
  lines.each {|gene| genes_list << gene.chomp}
  
return genes_list
end


#### OBTAIN THE INTERACTIONS OF THE GENES FROM THE ORIGINAL LIST
genes_list = read_file_genes(file)
hash_interactions = Hash.new
genes_list.each {|gene| 
hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 1)} #each value of the hash is a Gene object that correspond with the genes from the original list
 

#### STORE THE NEW GENES (Level 2) IN ANOTHER LIST, TRANSFORM THEM INTO A GENE OBJECT AND LOOK FOR THEIR INTERACTIONS
genes_list2 = Array.new()
hash_interactions.values.each {|gene_obj| gene_obj.interactions.each {|gene| genes_list2 << gene}}


genes_list2.each {|gene| 
hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 2)} #each value of the hash is a Gene object that correspond with the genes from the original list

#### STORE THE NEW GENES (Level 3) IN ANOTHER LIST, TRANSFORM THEM INTO A GENE OBJECT AND LOOK FOR THEIR INTERACTIONS
genes_list3 = Array.new()
hash_interactions.values[168..-1].each {|gene_obj|
  if gene_obj.level == 2
    gene_obj.interactions.each {|gene| genes_list3 << gene}
  else
    next
  end}

genes_list3.each {|gene|
   if hash_interactions.has_key?(gene)
       next   
   else 
     hash_interactions["#{gene.downcase}"] = Gene.new(:agi_locus => "#{gene.downcase}",:level => 3)
   end} #each value of the hash is a Gene object that correspond with the genes from the original list



### DELETING GENES THAT DON'T HAVE INTERACTIONS:
 
hash_interactions_updated = Hash.new #new hash that will contain only those genes that have interactions
hash_interactions.values.each{|gene_obj|
  if gene_obj.interactions.length < 1
    next
  end
  if gene_obj.level > 1 and gene_obj.interactions.length < 2
    next
  end
  hash_interactions_updated[gene_obj.agi_locus] = gene_obj}


######### CREATING AND OBTAINING INTERACTIONS NETWORKS USING THE INTERACTION NETWORK CLASS:
list_of_networks = Array.new
hash_interactions_updated.values.each {|gene|                            
                            if list_of_networks.any?{|int_network| int_network.network.any?{|int| gene.agi_locus == int.agi_locus}} #para no hacer redes repetidas
                              next
                            else
                              network = InteractionNetwork.new(:gene_original => gene,:all_genes => hash_interactions_updated)                               
                              list_of_networks << network
                            end}

#### A FUNCTION THAT STORE IN AN ARRAY ONLY THE NETWORKS THAT INCLUDES AT LEAST TWO OF THE ORIGINAL GENES:
def remove_networks_no_functional(array_of_networks,new_list_networks)
  array_of_networks.each {|int_network_obj| int_network = int_network_obj.network 
    if int_network.count{|gene| gene.level.to_i == 1} < 2
       next
    else
       new_list_networks << int_network_obj
    end}
  return new_list_networks
end

list_networks_updated = remove_networks_no_functional(list_of_networks,Array.new) #applies the function to the list of total networks obtained
puts "KEGG ANNOTATIONS: "
list_networks_updated.each{|net| puts net.kegg_annotations
 puts
 puts
 }
puts "GO ANNOTATIONS: "
list_networks_updated.each{|net| puts net.go_annotations
 puts
 puts
 } 
 
    