cat(paste0("current dir: ", getwd(), "\n"))

#' @description open/create a qmd or Rmd file and switch to the document directory
myopen <- function(file) {
  if (!tools::file_ext(file) %in% c("Rmd", "qmd")) {
    stop(cat("\nWrong file type. Only open/create .Rmd or .qmd files\n"))
  }

  current_path = dirname(".")
  doc_path = dirname(file)

  if (!file.exists(file)) {
    if (!dir.exists(doc_path)) dir.create(doc_path)
    file.create(file)

    yaml_title <- tools::file_path_sans_ext(basename(file))
    yaml_date <- as.character(Sys.Date())
    yaml_author <- ""
    git_config <- git2r::config()
    if (!is.null(git_config$global$user.name)) {
      yaml_author <- git_config$global$user.name
    }
    header <- glue::glue("---\ntitle: \"{yaml_title}\"  \nauthor: \"{yaml_author}\"  \ndate: \"{yaml_date}\"  \ncategories: [\'analysis\']  \ncode-fold: true  \nexecute:  \n  warning: false  \n  messages: false  \neditor_options:  \n  chunk_output_type: console  \n---")
    boiler <- c("", "## Introduction", "", "Write something ...")
    writeLines(c(header, boiler), file)
  }
  rstudioapi::navigateToFile(file)

  if (current_path != doc_path) {
    setwd(doc_path)
    cat(paste0("Set working dir to: \n", doc_path, "\n"))
  }
}


