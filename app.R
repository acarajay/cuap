# -----------------------------
# Author: Almie Carajay
# Date Started: June 30, 2017
# -----------------------------
library(shiny)
library(ape)
library(seqinr)
library(hash)
library(dplyr)
library(plyr)
library(stringr)
library(tidyr)
library(boot)   #Bootstrap functions
#Cluster, factoextra, magrittr,NbClust package for generating heat maps
library(cluster)
library(factoextra)
library(magrittr)
library(NbClust)
library(gplots)


shinyApp(
  
    ui = shinyUI(htmlTemplate("www/index.html")),
  
    server = shinyServer(function(input, output) {
    
	    data1 <- data.frame(iris)
	    output$plot1 <- renderPlot(plot(data1[,c(2,3)]))

		radiodata <- reactive({
			radioval <- input$radioval
			if(is.null(radioval)){return()}
			radioval
		})


	    output$sumofdata <- renderDataTable({
	    	radioval <- input$radioval

	    	#radioval <- radiodata() 
	    	#print(radioval)

	    	if(is.null(radioval)){return()}
	    	gene.Name <- ''

	    	summarydata <- function(seqlist, numofcodonList, acc){  # prints table of sequences and the number of codons analyzed
	    		j = 1
	    		annotList=''
	    		nrow=''    		

		        for(j in 1:length(seqlist)){
		        	nrow <- append(nrow, j)
		        	if(acc==1){
						annotList[j] <- c(gene.Name[j])
		        	}else{
		           		annotList[j] <- c(getAnnot(seqlist[[j]]))
		        	}
		            
		        }
		        frow <- nrow[2:length(nrow)]  

	    		sumtable <- data.frame(frow, annotList,numofcodonList)  # table of codons with frequency per codon
	    		colnames(sumtable) <- c("#","Sequence(s)", "Number of Codons Analyzed")
				sumtable
			}  

	    	codonrowcount <- function(seqlist, acc){
	    		#for each sequence in seqlist extract each codon
	        	a = 1
	        	#print(length(seqlist))
	        	dnaletters=''
	        	numofcodonList=''
	        	allcodons=''
	        	allaminoacids=''
	        	allidentifier=''
	        	allaminonames=''
	        	preferedindex=''
	        	preferedcodon=''

	        	#print(acc)
	        	#print(attr(seqlist,"description"))
	        	#print(speciesAccessno[c])
	        	for(a in 1:length(seqlist)){   #index of each sequence in seqlist
	        		#print(seqlist[a])  # but sequences here is of string value
	          		#extract each character in a sequence of seqlist
	          		invalid=FALSE

	          		
	          
	          		if(acc==1){
	          			seqstring <- trimSpace(seqlist[[a]])   # trims (" Some text. " == "Some text.")
	          		}else{
	          			seqstring <- trimSpace(seqlist[a])   # trims (" Some text. " == "Some text.")
	          		}
	          		newseq <- str_replace_all(seqstring, "[ \r\n]" , "")   #replace all newline, tabs and whitespace with ""
	          		bases <- unlist(str_split(newseq, ""))  #[[1]] splits into individual letters/characters
	          		dnaletters <- toupper(bases)

	          		#extract each codon for each sequence in seqlist
			        #print(length(dnaletters))
			        b = 1
			        codon=''
			        aminoacid=''
			        id_amino=''
			        amino_name=''
			        codonlist <- hash("ATA"="I","ATC"="I","ATT"="I","ATG"="M",
			                            "ACA"="T","ACC"="T","ACG"="T","ACT"="T",
			                            "AAC"="N","AAT"="N","AAA"="K","AAG"="K",
			                            "AGC"="S","AGT"="S","AGA"="R","AGG"="R",
			                            "CTA"="L","CTC"="L","CTG"="L","CTT"="L",
			                            "CCA"="P","CCC"="P","CCG"="P","CCT"="P",
			                            "CAC"="H","CAT"="H","CAA"="Q","CAG"="Q",
			                            "CGA"="R","CGC"="R","CGG"="R","CGT"="R",
			                            "GTA"="V","GTC"="V","GTG"="V","GTT"="V",
			                            "GCA"="A","GCC"="A","GCG"="A","GCT"="A",
			                            "GAC"="D","GAT"="D","GAA"="E","GAG"="E",
			                            "GGA"="G","GGC"="G","GGG"="G","GGT"="G",
			                            "TCA"="S","TCC"="S","TCG"="S","TCT"="S",
			                            "TTC"="F","TTT"="F","TTA"="L","TTG"="L",
			                            "TAC"="Y","TAT"="Y","TAA"="_","TAG"="_",
			                            "TGC"="C","TGT"="C","TGA"="_","TGG"="W")

			        identifier <- hash("ATA"="Ile3","ATC"="Ile2","ATT"="Ile1","ATG"="Met",
	                          "ACA"="Thr3","ACC"="Thr2","ACG"="Thr4","ACT"="Thr1",
	                          "AAC"="Asn2","AAT"="Asn1","AAA"="Lys1","AAG"="Lys2",
	                          "AGC"="Ser6","AGT"="Ser5","AGA"="Arg5","AGG"="Arg6",
	                          "CTA"="Leu3","CTC"="Leu2","CTG"="Leu4","CTT"="Leu1",
	                          "CCA"="Pro3","CCC"="Pro2","CCG"="Pro4","CCT"="Pro1",
	                          "CAC"="His2","CAT"="His1","CAA"="Gln1","CAG"="Gln2",
	                          "CGA"="Arg3","CGC"="Arg2","CGG"="Arg4","CGT"="Arg1",
	                          "GTA"="Val3","GTC"="Val2","GTG"="Val4","GTT"="Val1",
	                          "GCA"="Ala3","GCC"="Ala2","GCG"="Ala4","GCT"="Ala1",
	                          "GAC"="Asp2","GAT"="Asp1","GAA"="Glu1","GAG"="Glu2",
	                          "GGA"="Gly3","GGC"="Gly2","GGG"="Gly4","GGT"="Gly1",
	                          "TCA"="Ser3","TCC"="Ser2","TCG"="Ser4","TCT"="Ser1",
	                          "TTC"="Phe2","TTT"="Phe1","TTA"="Leu5","TTG"="Leu6",
	                          "TAC"="Tyr2","TAT"="Tyr1","TAA"="STOP","TAG"="STOP",
	                          "TGC"="Cys2","TGT"="Cys1","TGA"="STOP","TGG"="Trp")

			        aminoNames <- hash("ATA"="Isoleucine","ATC"="Isoleucine","ATT"="Isoleucine","ATG"="Methionine",
	                          "ACA"="Threonine","ACC"="Threonine","ACG"="Threonine","ACT"="Threonine",
	                          "AAC"="Asparagine","AAT"="Asparagine","AAA"="Lysine","AAG"="Lysine",
	                          "AGC"="Serine","AGT"="Serine","AGA"="Arginine","AGG"="Arginine",
	                          "CTA"="Leucine","CTC"="Leucine","CTG"="Leucine","CTT"="Leucine",
	                          "CCA"="Proline","CCC"="Proline","CCG"="Proline","CCT"="Proline",
	                          "CAC"="Histidine","CAT"="Histidine","CAA"="Glutamine","CAG"="Glutamine",
	                          "CGA"="Arginine","CGC"="Arginine","CGG"="Arginine","CGT"="Arginine",
	                          "GTA"="Valine","GTC"="Valine","GTG"="Valine","GTT"="Valine",
	                          "GCA"="Alanine","GCC"="Alanine","GCG"="Alanine","GCT"="Alanine",
	                          "GAC"="Aspartic acid","GAT"="Aspartic acid","GAA"="Glutamic acid","GAG"="Glutamic acid",
	                          "GGA"="Glycine","GGC"="Glycine","GGG"="Glycine","GGT"="Glycine",
	                          "TCA"="Serine","TCC"="Serine","TCG"="Serine","TCT"="Serine",
	                          "TTC"="Phenylalanine","TTT"="Phenylalanine","TTA"="Leucine","TTG"="Leucine",
	                          "TAC"="Tyrosine","TAT"="Tyrosine","TAA"="STOP CODON","TAG"="STOP CODON",
	                          "TGC"="Cysteine","TGT"="Cysteine","TGA"="STOP CODON","TGG"="Tryptophan")

					
					for(b in 1:length(dnaletters)){
	          			#print(dnaletters[b])    # "A","T","G","C"
	          			if(dnaletters[b]=="A" || dnaletters[b]=="T" || dnaletters[b]=="G" || dnaletters[b]=="C"){
	          				if(b%%3==0){   # every 3rd character
		              			nucleotide <- c(dnaletters[(b-2):b])    #"A""T""G"
		              			y <- paste(nucleotide, collapse="")  #"ATG"
		              			codon <- append(codon, y)
		              			aminoacid <- append(aminoacid, codonlist[[y]])              #prints list of protein ("M","M""V","M""V""S")
		                		id_amino <- append(id_amino, identifier[[y]])      # prints list of identifier
		                		amino_name <- append(amino_name, aminoNames[[y]])   # prints list names of Amino acid
		                	}
		               	}else{
	          				print(c("Sequence contains invalid input:", dnaletters[b]))
	          				invalid = TRUE
	          				break
	          			}	            
	          		}
	          

	            	if(invalid==TRUE){
	          			break
	          		}else{        
		        	#solve for the codon count 
		          
		        		proteinseq <- paste(aminoacid, collapse="")
		        		if(codon[1]=="" || protein[1]==""){
		        			codon <- codon[2:length(codon)]
		        			aminoacid <- aminoacid[2:length(aminoacid)]
		        			id_amino <- id_amino[2:length(id_amino)]
		        			amino_name <- amino_name[2:length(amino_name)]
		        		}
		        		#print(id_amino)
		        		#tt <- data.frame( codon, id_amino, aminoacid)
		        		#View(tt)
				        tb <- data.frame( codon, aminoacid)  # table of codons with frequency per codon

				        #View(tb)
				        numofcodon <- length(codon)
				        codoncount <- count(tb, c('codon','aminoacid'))   #counts for the frequency
				        #print(codoncount)
				        #View(codoncount)
				        #print(numofcodon)
				        #----------write into csv file the frequency counts ------------
		        		write.table(codoncount, file="widefinaldata.csv", sep=",", row.names=FALSE)
		                numofcodonList <- append(numofcodonList, numofcodon)
			        }
	            	allcodons <- append(allcodons,codon)  #appends each codons per sequence in a given input
	          		allaminoacids <- append(allaminoacids,aminoacid)
	          		allidentifier <- append(allidentifier, id_amino)
	          		allaminonames <- append(allaminonames, amino_name)   
		        }
		        CODONS <- allcodons[2:length(allcodons)]
		        AMINO_ACID <- allaminoacids[2:length(allaminoacids)]
		        IDENTIFIER <- allidentifier[2:length(allidentifier)]
		        NAME <- allaminonames[2:length(allaminonames)]
		        #print(CODONS)   #prints all codons in a given input except null values
		        #print(AMINO_ACID)
		        cut <- data.frame(CODONS, IDENTIFIER, AMINO_ACID, NAME)
		        newcut <- count(cut,c('CODONS', 'IDENTIFIER', 'AMINO_ACID','NAME'))
		        #View(newcut)
		        #print(length(newcut$CODONS))
		        
		        totalcount <- count(AMINO_ACID)
		        #print(totalcount)
		        #print(totalcount[1,1])
		        #print(totalcount[1,2])
		        #print(CODONS)
		        #print(length(newcut[,3]))

		        n = 1
		        m = 1
		        TOTALCOUNT = ''
		        for(n in 1:length(newcut[,3])){
		        	#print(newcut[n,3])
		        	for(m in 1:length(totalcount[,1])){
		        		#print(totalcount[m,1])
		        		if(newcut[n,3]==totalcount[m,1]){
		        			TOTALCOUNT <- append(TOTALCOUNT, totalcount[m,2])
		        		}
		        	}
		        }
		        newcut$TOTAL <- TOTALCOUNT[2:length(TOTALCOUNT)]
		        newcut <- transform(newcut, X = freq/as.numeric(TOTAL))
		        #print(newcut$freq)
		        #print(newcut$TOTAL)
		        #print(data.matrix(newcut$TOTAL))
				View(newcut)

		        #Maximum Likelihood estimation of parameters of the MLE based on the observed data
				mle<-function(x) {
				  experiments <- rmultinom(n=1, size=1000, prob=x)
				  df=data.frame(experiments)/1000
				  mle<-df[,1]
				  se<-sqrt((mle*(1-mle))/1000)
				  return(cbind(mle,se))
				}
				mtable <- data.frame(mle(newcut$X))
				View(mtable)
				codontable <- cbind(newcut, mtable)
				View(codontable)
				#codontable <- merge(newcut,mtable, by.x="X", by.x="y")

				#print(length(totalcount[,1]))
				f = 1
				g = 1
				
				LSElist=''
				
				#print(length(codontable[,9]))
				for(f in 1:length(totalcount[,1])){
					count = 0
					SElist =''
					indexlist=''
					for(g in 1:length(codontable[,3])){						
						if(totalcount[f,1]==codontable[g,3]){
							count = count + 1
							#print(c("Count:",count))
							SElist <- append(SElist, codontable[g,9])
							indexlist <- append(indexlist, g)
							#print("*******")
							#print(g)
							#print("*******")
							#print(SElist)
						}
					}
					#print("====================================")
					#print(count)

					newSElist <- SElist[2:length(SElist)]
					newindexlist <- indexlist[2:length(indexlist)]
					#print(newSElist)
					#print(newindexlist) 
					#print(totalcount[f,1])
					lse <- min(SElist[2:length(SElist)])

					h = 1					
					for(h in 1:length(newindexlist)){
						#print(newSElist[h])
						if(lse == newSElist[h]){
							#print("***************")
							e <- newindexlist[h]
							#print(newindexlist[h])
							preferedcodon <- append(preferedcodon, codontable[e,1])
							preferedindex <- append(preferedindex,e)
							#print(codontable[e,1])
							#print("***************")

							#newcodonlist <- append(newcodonlist, codontable[h,1])
						}
					}
					#print("====================================")
					#print(min(SElist[2:length(SElist)]))
					#print("====================================")
					LSElist <- append(LSElist, lse)
					
				}
				#print(preferedindex[2:length(preferedindex)])
				print("prefered index are: -----------")
				index <- as.numeric(preferedindex[2:length(preferedindex)])
				sortIndex <- sort(index, decreasing = FALSE)
				print(sortIndex)
				#PREFERRED_CODON <- preferedcodon[2:length(preferedcodon)]
				print("----Preferred codons:---------")
				#print(PREFERRED_CODON)
				p = 1
				q = 1
				for(p in 1:length(codontable[,1])){
					for(q in 1:length(sortIndex)){
						if(p == sortIndex[q]){
							preferedcodon <- append(preferedcodon, codontable[p,1])
							#print(codontable[p,1])
						}
					}

				}
				print(codontable[[1]][1])
				#print(CODONS)
				print(preferedcodon)

				LSE <- LSElist[2:length(LSElist)]
				#print(LSElist[2:length(LSElist)])
				AMINO <- totalcount[,1]
				LSEtable <- data.frame(AMINO,LSE)
				View(LSEtable)



		        #----------write into csv file the frequency counts ------------
		        #tdata <- reactive({write.table(newcut, file="tempdata.csv", sep=",", row.names=FALSE)})    
		        #wds<- data.frame(spread(newcut, CODONS, freq, fill=0, convert = TRUE))     #
				#View(wds)

				#ds1 <- data.matrix(wds)
				#ala1<-sum(ala$Ala1)/sum(ala[,1:4])
				#ala <- ds1[,16:19]
				#print(ala)
				#ala1 <- rowSums(ala[,1:4])
				#print(ala1)
				
				
				#Computing for the relative frequencies for alanine
				#ala$Rfala1<-ala$Ala1/rowSums(ala[,1:4])
				#print(length(AMINO_ACID))
				#function for generating heatmaps for the codon usage
				heatmap<-function(data) {
				    my_data<-data
					res.dist<-get_dist(my_data,stand=TRUE,method="spearman")
					#Visualize the dynamics of the distances
					fviz_dist(res.dist,gradient=list(low="#00AFBB",mid="white",high="#FC4E07"))
					#High "Orange"; LOw "Blue
				}
				#hm = heatmap(ds1)
				#print(hm)
				#aminotb <- data.frame(AMINO_ACID)
				#rowcount <- count(aminotb, c('AMINO_ACID'))    # sum of row count per amino acid
				#rowtb <- count(aminotb,c('AMINO_ACID','CODONS'))
				#View(rowcount)
				#sumeachrow <- apply(wds, 1, sum)
				#print(sumeachrow)
		        st <- summarydata(seqlist, numofcodonList[2:length(numofcodonList)], acc)
	    	}

	

	    	if(radioval=="one"){                 #enter sequence(s)
	    		genome <- input$sequence
	    		#----------------- getting sequence(s) from user input ------------------
	    		if(is.null(genome) || genome==""){return()}else{
			        #write input sequence into genome.fasta
			        writeGenome <- write.fasta(sequences = genome, names = names(genome), nbchar = 60, file.out = "genome.fasta")
			        #read the input sequence from a genome.fasta, sequence as DNA, string, uppercase, sequence attributes should be set, only sequences as returned , the '>' at the beginning of the description lines is removed in the annotations of the sequences
			        genefile <- read.fasta(file = "genome.fasta", seqtype = c("DNA"), as.string = TRUE, forceDNAtolower = FALSE, set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = TRUE)
	    		}
	    		#------------------------------------------------------------------------
	    		#create a new list containing only the sequences
		        seqlist=''
		        annotList=''
		        i = 2
		        j = 1
		        for(i in 2:length(genefile)){
		          seqlist[j] <- c(genefile[i])
		          annotList[j] <- c(getAnnot(genefile[[i]]))
		          j = j+1
		        }
		        acc <- 0
		        nc <- codonrowcount(seqlist, acc)

	    	}else if(radioval=="two"){
	    		access.nb <- input$accessno

	    		#----------------- getting sequence(s) from user input ------------------
	    		if(is.null(access.nb) || access.nb==""){return()}else{
			      	genefile <- read.GenBank(access.nb)
			        seqlist <- read.GenBank(access.nb, species.names = TRUE,
	                   gene.names = FALSE, as.character = TRUE)
			      #print(genefile)
                  #speciesName <- attr(genefile, "species") ## get the species names
			      #speciesAccessno <- names(genefile)
                  ## build a matrix with the species names and the accession numbers:
			      #gen <- cbind(attr(genefile, "species"), names(genefile))		      
			      gene.Name <- attr(genefile, "description")  ## the description from each FASTA sequence:
			      #print(gen)
	    		}
	      		#-----------------------------------------------------------------------
	      		acc <- 1
	      		st <- codonrowcount(seqlist, acc)

	    	}else if(radioval=="three"){
			    # input$file1 will be NULL initially. After the user selects
			    # and uploads a file, head of that data file by default,
			    # or all rows if selected, will be shown.
			    inputfile <- req(input$inputfile)

	    		if(is.null(inputfile)){return()}else{
			      # read the content of input file
			      genefile <- read.fasta(file = input$inputfile$datapath, seqtype = c("DNA"), as.string = TRUE, forceDNAtolower = FALSE, set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = TRUE)
			      #print(dnafile)
			      #print(length(dnafile))
			      acc <- 0
			      nc <- codonrowcount(genefile, acc)
	      		}

	    	}else{
	    		print(radioval)
	    		print("No sequence input!")
	    	}


	    })

		output$codonusagetable <- renderDataTable({
			sdata<-read.csv("tempdata.csv",header=TRUE, sep=",", fill=TRUE)
			#print(sdata)



		})	    

		
		
    })
  
)
