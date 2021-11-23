################
## ALIGNMENTS ##
##      &     ##
##    TREES   ##
################ 

Seq_HoxD <- readDNAStringSet ("02_Row_Data/HoxD_Sequences.fasta") # Load sequences 
print (Seq_HoxD) # Print

class (Seq_HoxD) # Check the class of the object 

## Alignment with Clustal W, Clustal O & Muscle 
Alig_HoxD_Clustal_W <- msa (Seq_HoxD, method = "ClustalW") # Make the alignment with Clustal W
Alig_HoxD_Clustal_O <- msa (Seq_HoxD, method = "ClustalOmega") # Make the alignment with Clustal o
Alig_HoxD_Muscle <- msa (Seq_HoxD, method = "Muscle") # Make the alignment with Muscle

print (Alig_HoxD_Clustal_W) # Print 
print (Alig_HoxD_Clustal_O) # Print 
print (Alig_HoxD_Muscle) # Print

T_Clustal_W <- msaConvert (Alig_HoxD_Clustal_W, type = "seqinr::alignment") # Convert the data for the sequence analysis package "seqinr" 
T_Clustal_O <- msaConvert (Alig_HoxD_Clustal_O, type = "seqinr::alignment") # Convert the data for the sequence analysis package "seqinr"
T_Muscle <- msaConvert (Alig_HoxD_Muscle, type = "seqinr::alignment") # Convert the data for the sequence analysis package "seqinr"

Dist_Clustal_W <- dist.alignment (T_Clustal_W, "identity") # Create the matrix of distances from the alignment, an identity matrix
as.matrix (Dist_Clustal_W) # Save it as matrix 

Dist_Clustal_O <- dist.alignment (T_Clustal_O, "identity") # Create the matrix of distances from the alignment, an identity matrix
as.matrix (Dist_Clustal_O) # Save it as matrix 

Dist_Muscle <- dist.alignment (T_Muscle, "identity") # Create the matrix of distances from the alignment, an identity matrix
as.matrix (Dist_Muscle) # Save it as matrix 


## Phylogenetic trees
Tree_Clustal_W <- njs (Dist_Clustal_W) # Make the tree with neighbor joining
print (Tree_Clustal_W) # Print

Tree_Clustal_O <- njs (Dist_Clustal_O) # Make the tree with neighbor joining
print (Tree_Clustal_O) # Print

Tree_Muscle <- nj (Dist_Muscle) # Make the tree with neighbor joining
print (Tree_Muscle) # Print


# Simple trees
Simple_Tree_Clustal_W <- ggtree (Tree_Clustal_W) # Make the information a tree
pdf (file = "03_Output/Plots/Simple_Tree_Clustal_W.pdf", width = 10, height = 10) # Save it as PDF with specific measures 
plot (Simple_Tree_Clustal_W) # Run the plot
dev.off () # Close the output

Simple_Tree_Clustal_O <- ggtree (Tree_Clustal_O) # Make the information a tree
pdf (file = "03_Output/Plots/Simple_Tree_Clustal_O.pdf", width = 10, height = 10) # Save it as PDF with specific measures 
plot (Simple_Tree_Clustal_O) # Run the plot
dev.off () # Close the output

Simple_Tree_Muscle <- ggtree (Tree_Muscle) # Make the information a tree
pdf (file = "03_Output/Plots/Simple_Tree_Muscle.pdf", width = 10, height = 10) # Save it as PDF with specific measures 
plot (Simple_Tree_Muscle) # Run the plot
dev.off () # Close the output


## Not so simple trees
Color_Tree_Clustal_W <- ggtree (Tree_Clustal_W, aes (color = branch.length)) + # Make the information a tree, colored by the branch lengths
  xlim (0, 1) + geom_tiplab (size = 4, color = "darkblue") +  # Determine the limits and add the tip label layer, with blue color
  geom_label2 (aes (subset = !isTip, label = node), size = 4, color = "darkblue", alpha = 0.5) + # Support the aesthetic: labels of nodes, color white and transparency of 0.5
  theme_tree2 ("white") + theme (legend.position = "bottom") # Use the tree2 theme (background) = black, and place the legend on the bottom
pdf ("03_Output/Plots/Color_Tree_Clustal_W.pdf", width = 10, height = 10) # Save it as PDF with specific measures 
plot (Color_Tree_Clustal_W) # Run the plot
dev.off () # Close the output

Color_Tree_Clustal_O <- ggtree (Tree_Clustal_O, aes (color = branch.length)) + # Make the information a tree, colored by the branch lengths
  xlim (0, 1) + geom_tiplab (size = 4, color = "darkblue") +  # Determine the limits and add the tip label layer, with blue color
  geom_label2 (aes (subset = !isTip, label = node), size = 4, color = "darkblue", alpha = 0.5) + # Support the aesthetic: labels of nodes, color white and transparency of 0.5
  theme_tree2 ("white") + theme (legend.position = "bottom") # Use the tree2 theme (background) = black, and place the legend on the bottom
pdf ("03_Output/Plots/Color_Tree_Clustal_O.pdf", width = 10, height = 10) # Save it as PDF with specific measures 
plot (Color_Tree_Clustal_O) # Run the plot
dev.off () # Close the output

Color_Tree_Muscle <- ggtree (Tree_Muscle, aes (color = branch.length)) + # Make the information a tree, colored by the branch lengths
  xlim (0, 0.25) + geom_tiplab (size = 4, color = "darkblue") +  # Determine the limits and add the tip label layer, with blue color
  geom_label2 (aes (subset = !isTip, label = node), size = 4, color = "darkblue", alpha = 0.5) + # Support the aesthetic: labels of nodes, color white and transparency of 0.5
  theme_tree2 ("white") + theme (legend.position = "bottom") # Use the tree2 theme (background) = black, and place the legend on the bottom
pdf ("03_Output/Plots/Color_Tree_Muscle.pdf", width = 10, height = 10) # Save it as PDF with specific measures 
plot (Color_Tree_Muscle) # Run the plot
dev.off () # Close the output
