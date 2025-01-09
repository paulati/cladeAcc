library(glue)

load_config <- function() {

    config_file_path <- system.file('config.yml', package = 'cladeAcc',
                                    mustWork = TRUE)
    config_default <- config::get(file = config_file_path, use_parent = FALSE,
                          config = 'default')

    Sys.setenv(R_CONFIG_ACTIVE = config_default$env)

    config <- config::get(file = config_file_path, use_parent = FALSE,
                                  config = Sys.getenv('R_CONFIG_ACTIVE'))
    return(config)
}

load_aws_config <- function() {

    config_file_path <- system.file('s3.yml', package = 'cladeAcc',
                                    mustWork = FALSE)
    config_default <- config::get(file = config_file_path, use_parent = FALSE,
                                  config = 'default')

    return(config_default)
}


# https://rstudio.github.io/config/articles/introduction.html
# https://docs.ansible.com/ansible/latest/reference_appendices/YAMLSyntax.html

# TODO: ver si esta es llamada desde algun lado
get_source_urls <- function(config, alignment_id) {

    files_url <- lapply(config$data_preparation$source$chrs, function(chr) {
        x <- glue::glue(config$data_preparation$source$file_name_pattern)
        file.path(config$data_preparation$source$base_url, x)
    })

    return(files_url)

}

# get_source_md5_url <- function(config, alignment_id) {
#
#     result <- file.path(config$data_preparation$source$base_url,
#                         config$data_preparation$source$md5_file_name)
#     return(result)
# }

# config <- load_config()
