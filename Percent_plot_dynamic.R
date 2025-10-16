library(ggplot2)
library(jsonlite)

make_plot = function(json_file) {
  params_json = fromJSON(json_file)
  if (params_json$input_table$ext == "tabular") {
    sep = "\t"
  } else {
    sep = ","
  }
  df = read.delim(params_json$input_table$source_path, header = TRUE, sep = sep)
  selected_col = params_json$selected_column
  samples_to_include = params_json$selected_value

  print(samples_to_include)
  
  #df_sub = df[df$orig_ident%in%samples_to_include, ]
  df_sub = df[df[[selected_col]] %in% samples_to_include, ]

  sumDf = data.frame(table(df_sub[, selected_col], df_sub$Orig_sample))
  colnames(sumDf)[1:2] = c("Cluster", "Group")
  
  groups = unique(sumDf$Group)
  
  comb_df = data.frame()
  for (cur_group in groups) {
    cur_df = sumDf[sumDf$Group == cur_group,]
    cur_df$Percent_Cluster = (cur_df$Freq / sum(cur_df$Freq)) * 100
    comb_df = rbind(comb_df, cur_df)
  }
  
  
  cur_plot = ggplot(comb_df, aes(fill=Cluster, y=Percent_Cluster, x=Group)) + 
    geom_bar(position="stack", stat="identity")+
    theme_classic() +
    theme(text = element_text(size = 4),
          axis.text.x = element_text(angle = 14, hjust = 1))
  
  #cur_plot
  
  ggsave("Percent_plot.png", cur_plot,  width = 2, height = 3, dpi = 300, units = 'in') 
  
}


if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  make_plot(args[1])
}

