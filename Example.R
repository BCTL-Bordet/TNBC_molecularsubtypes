## EXAMPLE

# 1. Load  gene signatures  
load("~/lehmann.RData")

# 2. Load  gene expression data  
load("~/normalized_counts.RData") ## The data was previosly normalized usin DESEQ2 vst. 

# 3. Calculate TNBC subtypes as chategories as defined by Bareche.
TNBCsubtype<- TNBCclassif(as.matrix(normalized_counts), version = "bareche", sig = sig, coef=F)
head(TNBCsubtype)
# Sample1 
# "immunomodulatory" 
# Sample2 
# "mesenchymal" 
# Sample3
# "basal_like_1" 
# Sample4
# "basal_like_1" 
# Sample5 
# "mesenchymal_stem_like" 
# Sample6 
# "mesenchymal_stem_like" 

# 4. Calculate TNBC subtypes with coefficients as defined by Bareche.
TNBCsubtype<- TNBCclassif(as.matrix(normalized_counts), version = "bareche", sig = sig, coef=T)
#         basal_like_1 immunomodulatory  mesenchymal mesenchymal_stem_like  luminal_ar
# Sample1     0.1575568       0.22035503 -0.006280037            0.06499573 -0.11533833
# Sample2    0.2491158      -0.26724748  0.362004149           -0.07653982 -0.12656296
# Sample3    0.4037095       0.03274814  0.033464367           -0.69405901 -0.34925315
# Sample4    0.1315930      -0.07942436  0.088592157           -0.31747438 -0.23614606
# Sample5   -0.4438737      -0.06607181  0.118596339            0.50761600  0.16024605
# Sample6   -0.1344419      -0.02239535 -0.017298300            0.10796467  0.08315154


# 5. Calculate TNBC subtypes as chategories as defined by Lehmann.
TNBCsubtype_leh<- TNBCclassif(as.matrix(normalized_counts), version = "lehmann", sig = sig, coef=F)
# Sample1 
# "immunomodulatory" 
# Sample2 
# "mesenchymal" 
# Sample3 
# "basal_like_1" 
# Sample4 
# "basal_like_1" 
# Sample5 
# "mesenchymal_stem_like" 
# Sample6 
# "basal_like_2" 

# 4. Calculate TNBC subtypes with coefficients as defined by Lehmann.
TNBCsubtype_leh<- TNBCclassif(as.matrix(normalized_counts), version = "lehmann", sig = sig, coef=T)

#           basal_like_1 basal_like_2 immunomodulatory  mesenchymal mesenchymal_stem_like  luminal_ar
# Sample1     0.1575568   0.08548929       0.22035503 -0.006280037            0.06499573 -0.11533833
# Sample2    0.2491158  -0.32440229      -0.26724748  0.362004149           -0.07653982 -0.12656296
# Sample3    0.4037095  -0.06984233       0.03274814  0.033464367           -0.69405901 -0.34925315
# Sample4    0.1315930  -0.21640421      -0.07942436  0.088592157           -0.31747438 -0.23614606
# Sample5   -0.4438737   0.43973493      -0.06607181  0.118596339            0.50761600  0.16024605
# Sample6   -0.1344419   0.17552633      -0.02239535 -0.017298300            0.10796467  0.08315154

