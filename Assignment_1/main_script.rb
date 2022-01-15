require "/home/osboxes/BioProCha_Assignments_LSG/Assignment_1/database_object.rb"
require "date"
#save the files name in different variables. The order has to be the same as in the commands we put in the shell
seed_stock_file,gene_file,hybrid_cross_file, seed_stock_modified = ARGV

#create a new object with the Database Class that will contain the three files' values
database = Database.new(:gene_file=>gene_file,:seed_stock_file => seed_stock_file,:hybrid_cross_file => hybrid_cross_file)

##### 1. PLANT 7 GRAMS OF SEEDS FROM EACH RECORD IN SEED STOCK AND UPDATE THE INFORMATION IN A NEW FILE:
##applying the modify_stock method to each seed object from the hash of seeds objects of the database object:

database.hash_seed_stock.values.each {|seed_object| seed_object.modify_stock(7)} #plants 7 grams of seeds

##modify the date of last planted by using the Class date of Ruby and the method Date.today:
database.hash_seed_stock.values.each {|seed_object| seed_object.last_planted = Date.today}

##creating the new file with the information of the seed stock database updated:
new_file = File.open("report.txt","w") #variable that will allows us to write the new file
header = IO.readlines(seed_stock_file)[0].chomp.split("\t") #obtains the header of the file (remove the "\n" at the end of the header with chomp)
header.each {|column| new_file << column + "\t"} #adds the columns names in the new file

#adds the new values of the columns using a loop: for each object, converts its properties' values into string and concatenates them. 
database.hash_seed_stock.values.each {|seed_object| new_file << "\n" + seed_object.seed_stock.to_s + "\t" + seed_object.mutant_gene_ID.Gene_ID.to_s +
   "\t" + seed_object.last_planted.to_s + "\t" + seed_object.storage.to_s + "\t" + seed_object.grams_remaining.to_s}
new_file.close #close the file

#### 2. WHICH GENES ARE GENETICALLY-LINKED? : CHI-SQUARED TEST

#applies the method chi_square to the cross objects of the hash cross, which is a property of the database object
database.hash_cross.values.each {|cross_object| cross_object.chi_square}
