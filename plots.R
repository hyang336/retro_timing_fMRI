library(rio)

data_yaam=import("C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\comp_table_yaam.csv")
data_haam=import("C:\\Users\\haozi\\OneDrive\\Desktop\\Postdoc\\Wagner\\retro_timing_fMRI_data\\comp_table_haam.csv")

unique(data_haam$sub)
unique(data_yaam$sub)

hist(data_yaam$rec_error[data_yaam$smoothed=='no'],breaks=50,xlim=c(-2,2))

hist(data_yaam$rec_error[data_yaam$smoothed=='yes'],breaks=50,xlim=c(-2,2))

hist(data_haam$rec_error[data_haam$smoothed=='no'],breaks=50,xlim=c(-2,2))

hist(data_haam$rec_error[data_haam$smoothed=='yes'],breaks=50,xlim=c(-2,2))
