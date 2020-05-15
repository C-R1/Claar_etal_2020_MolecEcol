# Load necessary data
load("analyses/16S_ITS2/KI_Sym_Microbes_16S_ITS2.RData")

# Calculate number of samples with each genus (prev. clade)
sum(sample_data(phyASV.mic)$clade_A=='Y') #14 ***********************
sum(sample_data(phyASV.mic)$clade_A=='N') #212

sum(sample_data(phyASV.mic)$clade_B=='Y') #2
sum(sample_data(phyASV.mic)$clade_B=='N') #224

sum(sample_data(phyASV.mic)$clade_D=='Y') #145 *************************
sum(sample_data(phyASV.mic)$clade_D=='N') #81

sum(sample_data(phyASV.mic)$clade_F=='Y') #2
sum(sample_data(phyASV.mic)$clade_F=='N') #224

sum(sample_data(phyASV.mic)$clade_G=='Y') #19 *********************
sum(sample_data(phyASV.mic)$clade_G=='N') #207


sum(sample_data(phyASV.mic.FaviaM)$clade_A=='Y') #0
sum(sample_data(phyASV.mic.FaviaM)$clade_A=='N') #8

sum(sample_data(phyASV.mic.FaviaM)$clade_B=='Y') #0
sum(sample_data(phyASV.mic.FaviaM)$clade_B=='N') #8

sum(sample_data(phyASV.mic.FaviaM)$clade_D=='Y') #4
sum(sample_data(phyASV.mic.FaviaM)$clade_D=='N') #4

sum(sample_data(phyASV.mic.FaviaM)$clade_F=='Y') #0
sum(sample_data(phyASV.mic.FaviaM)$clade_F=='N') #8

sum(sample_data(phyASV.mic.FaviaM)$clade_G=='Y') #0
sum(sample_data(phyASV.mic.FaviaM)$clade_G=='N') #8



sum(sample_data(phyASV.mic.Favites)$clade_A=='Y') #0
sum(sample_data(phyASV.mic.Favites)$clade_A=='N') #23

sum(sample_data(phyASV.mic.Favites)$clade_B=='Y') #0
sum(sample_data(phyASV.mic.Favites)$clade_B=='N') #23

sum(sample_data(phyASV.mic.Favites)$clade_D=='Y') #23
sum(sample_data(phyASV.mic.Favites)$clade_D=='N') #0

sum(sample_data(phyASV.mic.Favites)$clade_F=='Y') #0
sum(sample_data(phyASV.mic.Favites)$clade_F=='N') #23

sum(sample_data(phyASV.mic.Favites)$clade_G=='Y') #0
sum(sample_data(phyASV.mic.Favites)$clade_G=='N') #23

sum(sample_data(phyASV.mic.Hydnophora)$clade_A=='Y') #0
sum(sample_data(phyASV.mic.Hydnophora)$clade_A=='N') #38

sum(sample_data(phyASV.mic.Hydnophora)$clade_B=='Y') #0
sum(sample_data(phyASV.mic.Hydnophora)$clade_B=='N') #38

sum(sample_data(phyASV.mic.Hydnophora)$clade_D=='Y') #12 ****************
sum(sample_data(phyASV.mic.Hydnophora)$clade_D=='N') #26

sum(sample_data(phyASV.mic.Hydnophora)$clade_F=='Y') #0
sum(sample_data(phyASV.mic.Hydnophora)$clade_F=='N') #38

sum(sample_data(phyASV.mic.Hydnophora)$clade_G=='Y') #2
sum(sample_data(phyASV.mic.Hydnophora)$clade_G=='N') #36


sum(sample_data(phyASV.mic.Montipora)$clade_A=='Y') #1
sum(sample_data(phyASV.mic.Montipora)$clade_A=='N') #65

sum(sample_data(phyASV.mic.Montipora)$clade_B=='Y') #2
sum(sample_data(phyASV.mic.Montipora)$clade_B=='N') #64

sum(sample_data(phyASV.mic.Montipora)$clade_D=='Y') #42 ********************
sum(sample_data(phyASV.mic.Montipora)$clade_D=='N') #24

sum(sample_data(phyASV.mic.Montipora)$clade_F=='Y') #1
sum(sample_data(phyASV.mic.Montipora)$clade_F=='N') #65

sum(sample_data(phyASV.mic.Montipora)$clade_G=='Y') #11 **********************
sum(sample_data(phyASV.mic.Montipora)$clade_G=='N') #55


sum(sample_data(phyASV.mic.Platygyra)$clade_A=='Y') #10 *******************
sum(sample_data(phyASV.mic.Platygyra)$clade_A=='N') #22

sum(sample_data(phyASV.mic.Platygyra)$clade_B=='Y') #0
sum(sample_data(phyASV.mic.Platygyra)$clade_B=='N') #32

sum(sample_data(phyASV.mic.Platygyra)$clade_D=='Y') #32
sum(sample_data(phyASV.mic.Platygyra)$clade_D=='N') #0

sum(sample_data(phyASV.mic.Platygyra)$clade_F=='Y') #0
sum(sample_data(phyASV.mic.Platygyra)$clade_F=='N') #32

sum(sample_data(phyASV.mic.Platygyra)$clade_G=='Y') #1
sum(sample_data(phyASV.mic.Platygyra)$clade_G=='N') #31


sum(sample_data(phyASV.mic.Pocillopora)$clade_A=='Y') #0
sum(sample_data(phyASV.mic.Pocillopora)$clade_A=='N') #15

sum(sample_data(phyASV.mic.Pocillopora)$clade_B=='Y') #0
sum(sample_data(phyASV.mic.Pocillopora)$clade_B=='N') #15

sum(sample_data(phyASV.mic.Pocillopora)$clade_D=='Y') #8 ***********************
sum(sample_data(phyASV.mic.Pocillopora)$clade_D=='N') #7

sum(sample_data(phyASV.mic.Pocillopora)$clade_F=='Y') #0
sum(sample_data(phyASV.mic.Pocillopora)$clade_F=='N') #15

sum(sample_data(phyASV.mic.Pocillopora)$clade_G=='Y') #1
sum(sample_data(phyASV.mic.Pocillopora)$clade_G=='N') #14


sum(sample_data(phyASV.mic.Porites)$clade_A=='Y') #3
sum(sample_data(phyASV.mic.Porites)$clade_A=='N') #41

sum(sample_data(phyASV.mic.Porites)$clade_B=='Y') #0
sum(sample_data(phyASV.mic.Porites)$clade_B=='N') #44

sum(sample_data(phyASV.mic.Porites)$clade_D=='Y') #24 ******************
sum(sample_data(phyASV.mic.Porites)$clade_D=='N') #20

sum(sample_data(phyASV.mic.Porites)$clade_F=='Y') #1
sum(sample_data(phyASV.mic.Porites)$clade_F=='N') #43

sum(sample_data(phyASV.mic.Porites)$clade_G=='Y') #4
sum(sample_data(phyASV.mic.Porites)$clade_G=='N') #40
