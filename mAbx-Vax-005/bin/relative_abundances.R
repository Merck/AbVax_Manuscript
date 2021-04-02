# Copyright © 2021 Merck Sharp & Dohme Corp. a subsidiary of Merck & Co., Inc., Kenilworth, NJ, USA.
get_relative_count <- function(input){

  # Look at what these add up to
  tcount <- input %>%
  	group_by(names) %>%
  	summarize(total = sum(Count)) %>%
  	as.data.frame()

  inmerge <- inner_join(input, tcount, by = "names")
  inmerge$relabund<- 100 * inmerge$Count / inmerge$total

  # Confirm that they add to 100
  inmerge %>%
  	group_by(names) %>%
  	summarize(totalpercent = sum(relabund)) %>%
  	as.data.frame()

return(inmerge)
}
