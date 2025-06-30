Quantification using FeatureCount

#############QUANTIFICATION_PARAMETERS#########
# using infer_exp to determine strandedness. we use featureCounts from the subread package to quantify gene

# -T -> number of processors
# -a -> provide annotation file here
# -o -> provide outprefix
# -t -> provide feature type to quantify for
# -g -> attribute to assign the quantification for
# -s -> specify strandedness here, 0 = unstranded, 1 = stranded and 2 = reverse
# -p -> specify is data is paired end
# -B -> quanity if both ends of the pairs are mapped
