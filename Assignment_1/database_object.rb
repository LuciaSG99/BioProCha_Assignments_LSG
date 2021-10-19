class Gene #Object gene that allows us to modify the gene bank values 
   #Properties: (same as the columns of the file)
  attr_accessor :Gene_ID
  attr_accessor :Gene_name
  attr_accessor :mutant_phenotype
  
  def initialize (params = {}) #the values of the properties that the new object will have at the beginning 
    @Gene_ID = params.fetch(:Gene_ID,nil) 
    @Gene_name = params.fetch(:Gene_name,nil)
    @mutant_phenotype = params.fetch(:mutant_phenotype,nil)
  end
  
end

class Hybrid_Cross #Object Hybrid Cross that allows us to modify the cross bank values 
   #Properties: (same as the columns of the file)
  attr_accessor :Parent1
  attr_accessor :Parent2
  attr_accessor :F2_Wild
  attr_accessor :F2_P1
  attr_accessor :F2_P2
  attr_accessor :F2_P1P2
  
  def initialize (params = {}) #the values of the properties that the new object will have 
    @Parent1 = params.fetch(:Parent1,nil)
    @Parent2 = params.fetch(:Parent2,nil)
    @F2_Wild = params.fetch(:F2_Wild,nil)
    @F2_P1 = params.fetch(:F2_P1,nil)
    @F2_P2 = params.fetch(:F2_P2,nil)
    @F2_P1P2 = params.fetch(:F2_P1P2,nil)
  end
  
  def chi_square() #method that applies a chi-square test to the cross object
   total = @F2_Wild.to_f + @F2_P1.to_f + @F2_P2.to_f + @F2_P1P2.to_f #sum of all columns (properties). to_f : to transform them into float number
   
   #Calculate the expected observations:
   exp_f2w = 9.0/16.0*total.to_f 
   exp_f2p1 = 3.0/16.0*total.to_f
   exp_f2p2 = 3.0/16.0*total.to_f
   exp_f2p1p2 = 1.0/16.0*total.to_f
   
   #Calculate the chisquare value: (applying the formula)
   chi_square_value = ((@F2_Wild.to_f - exp_f2w)**2/exp_f2w)+((@F2_P1.to_f - exp_f2p1)**2/exp_f2p1)+((@F2_P2.to_f - exp_f2p2)**2/exp_f2p2)+((@F2_P1P2.to_f - exp_f2p1p2)**2/exp_f2p1p2) 
   if chi_square_value > 7.815 # chi-square value for a p < 0.05
    puts "Recording: #{@Parent1.mutant_gene_ID.Gene_name} is genetically linked to #{@Parent2.mutant_gene_ID.Gene_name} with chisquare score #{chi_square_value}"
   # The @Partent1 and @Parent2 have the  gene object inside (which is inside the stock object), which allow us to access to the property "Gene_name"
   end
   
  end
end 
 
class SeedStock #Object SeedStock allows us to modify the seed stock bank values 
  #Properties: (same as the columns of the file)
  attr_accessor :seed_stock
  attr_accessor :grams_remaining
  attr_accessor :mutant_gene_ID
  attr_accessor :last_planted
  attr_accessor :storage
  
  def initialize (params = {}) #the values of the properties that the new object will have 
    @grams_remaining = params.fetch(:grams_remaining,nil)
    @seed_stock = params.fetch(:seed_stock,nil)
    @mutant_gene_ID = params.fetch(:mutant_gene_ID,nil)
    @last_planted = params.fetch(:last_planted,nil)
    @storage = params.fetch(:storage,nil)
  end
  
  def modify_stock (n_seeds_planted) #method that allows us to change the number of seeds of our database
    @grams_remaining = @grams_remaining - n_seeds_planted 
    if @grams_remaining <= 0 #if there is no more seeds in the stock:
       puts "WARNING! we have run out of Seed Stock #{@seed_stock}"
       @grams_remaining = 0 #to force that the negative values are zero because we cannot have negative number of grams of seeds.
    end
    
    return @grams_remaining #return the new value of the property
    
  end
  
end


class Database # class that represents all the elements in our database
 attr_accessor :hash_seed_stock #hash that will contain all the SeedStock objects obtain from the file. 
 attr_accessor :hash_genes #hash that will contain all the Gene objects obtain from the file.
 attr_accessor :hash_cross #hash that will contain all the Hybrid_Cross objects obtain from the file.
 
def initialize (params = {})#the values of the properties that the new object database will have
  @hash_genes = get_genes(params.fetch(:gene_file)) #this property is the result of applying the function created below.
  #It's a hash where the values are genes objects and the keys are the genes ID
  @hash_seed_stock = get_seed_stock(params.fetch(:seed_stock_file))#same than before, but the values of the hash are SeedStock objects  
  @hash_cross = get_cross(params.fetch(:hybrid_cross_file))# the values of the hash are HYbrid_Cross objects
  
end
 
 #method that allows us to obtain the values of the table gene, store them as properties of a Gene object
 #and then all the objects are store in the hash_genes. It just needs as a argument the name of the file:
 
 def get_genes(gene_file) 
   lines = IO.readlines(gene_file) #array that contains the lines of the file
   hash_genes = Hash.new() #hash that will contain the Gene objects
   lines[1..-1].each {|line| gene_ID,gene_name, mutant_phenotype = line.split("\t") #iteraction where each value of a line is store in a variable
   #then it stores in the hash the new element created. The values of its properties will be the values of the previous variables.
   #The keys are the Gene ID (as a symbol)
   hash_genes[gene_ID.to_sym] = Gene.new(:Gene_ID=> gene_ID,:Gene_name =>gene_name,:mutant_phenotype =>mutant_phenotype)}
 
   return hash_genes #returns a hash where every element is a gene (as a object with its respective properties) 
 end
 
 #method that allows us to obtain the values of the table seed stock, store them as properties of a SeedStock object
 #and then all the objects are store in the hash_seeds. It just needs as a argument the name of the file:
 
 def get_seed_stock(seed_stock_file)
   lines = IO.readlines(seed_stock_file) #array that contains the lines of the file
   hash_seeds = Hash.new() #hash that will contain the SeedStock objects
   lines[1..-1].each {|line| seed_stock,mutant_gene_ID,last_planted,storage,grams_remaining = line.chomp.split("\t")#iteraction where each value of a line is store in a variable
   #then it stores in the hash the new element created. The values of its properties will be the values of the previous variables, with the exception of the
   #attribute mutant_gene_ID that will content the Gene object correspondent to that ID (because mutant_gene_ID and Gene_ID have the same value). 
   #The keys are the seed stock (as a symbol)
   hash_seeds[seed_stock.to_sym] = SeedStock.new(:seed_stock=> seed_stock,:mutant_gene_ID => @hash_genes[mutant_gene_ID.to_sym],:last_planted =>last_planted,:storage => storage,:grams_remaining => grams_remaining.to_i)}
  
   return hash_seeds #returns a hash where every element is a seed stock (as a object with its respective properties)
 end
 
 
 #method that allows us to obtain the values of the table Cross, store them as properties of a Hybrid_Cross object
 #and then all the objects are store in the hash_cross. It just needs as a argument the name of the file:
 
 def get_cross(hybrid_cross_file)
   lines = IO.readlines(hybrid_cross_file)#array that contains the lines of the file
   hash_cross = Hash.new() #hash that will contain the SeedStock objects
   lines[1..-1].each {|line| parent1,parent2,f2_wild,f2_P1,f2_P2,f2_P1P2 = line.split("\t")#iteraction where each value of a line is store in a variable
   #then it stores in the hash the new element created. The values of its properties will be the values of the previous variables, with the exception of 
   #attributes Parent1 and Parent2 that will content the Stock object correspondent to that seed stock (because this columns have as values seed stocks). 
   #The keys are --> parent1_parent2 (as a symbol)
   hash_cross[(parent1 + "_" + parent2).to_sym]= Hybrid_Cross.new(:Parent1=> @hash_seed_stock[parent1.to_sym],:Parent2 =>@hash_seed_stock[parent2.to_sym],:F2_Wild =>f2_wild,:F2_P1=> f2_P1,:F2_P2 => f2_P2,:F2_P1P2 => f2_P1P2.chomp)}
  
   return hash_cross  #returns a hash where every element is a cross (as a object with its respective properties)
 end
  
end
 
 


  





  
  
  
  
  
  
  
  
  
