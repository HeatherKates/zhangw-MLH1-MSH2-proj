read_HOMER <- function(file_path) {
  # Step 1: Read the raw text from the file
  lines <- readLines(file_path)
  
  # Step 2: Use regular expression to isolate the first column (Motif Name)
  motif_names <- sub("\t.*", "", lines)  # Extract text before the first tab character
  
  # Step 3: Remove the first column (Motif Name) and leave the rest of the line
  remaining_columns <- sub("^[^\t]*\t", "", lines)  # Remove the first column and the first tab
  
  # Step 4: Now split the remaining columns by tabs
  remaining_data_split <- strsplit(remaining_columns, "\t")
  
  # Step 5: Combine the first column (Motif Name) with the rest of the data
  remaining_data_df <- do.call(rbind, remaining_data_split)
  
  # Step 6: Combine the Motif Names with the rest of the columns
  col_names <- c("Motif Name", "Consensus", "P-value", "Log P-value", "q-value (Benjamini)",
                 "# of Target Sequences with Motif", "% of Target Sequences with Motif", 
                 "# of Background Sequences with Motif", "% of Background Sequences with Motif")
  
  df <- data.frame(motif_names, remaining_data_df, stringsAsFactors = FALSE)
  colnames(df) <- col_names
  
  # Step 7: Remove the '%' sign and convert the percentage columns to numeric
  df$`% of Target Sequences with Motif` <- as.numeric(gsub("%", "", df$`% of Target Sequences with Motif`))
  df$`% of Background Sequences with Motif` <- as.numeric(gsub("%", "", df$`% of Background Sequences with Motif`))
  
  # Step 8: Convert the numeric columns to the correct data types
  df[, 3:9] <- lapply(df[, 3:9], as.numeric)
  
  # Step 9: Remove the first line if necessary (you mentioned skipping the first row)
  df <- df[-1, ]
  
  # Return the cleaned data frame
  return(df)
}

