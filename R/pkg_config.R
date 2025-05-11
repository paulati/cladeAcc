
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
