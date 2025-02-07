library(aws.s3)
library(R.utils)

# https://cran.r-project.org/web/packages/heapsofpapers/vignettes/s3_setup.html


# bucketlist(key = "s3_upload", secret = "", add_region = "us-east-1")

# private functions

configure_aws_account <- function(access_key_id = NA, secret_acces_key = NA,
                                  default_region = NA)  {

    aws_config <- load_aws_config()
    if(is.na(access_key_id)) {
        access_key_id <- aws_config$aws_access_key_id
    }

    if(is.na(secret_acces_key)) {
        secret_acces_key <-aws_config$aws_secret_access_key
    }

    if(is.na(default_region)) {
        default_region <- aws_config$aws_default_region
    }

    Sys.setenv(
        "AWS_ACCESS_KEY_ID" = access_key_id,
        "AWS_SECRET_ACCESS_KEY" = secret_acces_key,
        "AWS_DEFAULT_REGION" = default_region
    )

}




# ------------------


download_from_s3 <- function(account_key, account_secret,
                             remote_base_folder_path, file_name_gz,
                             local_base_folder_path, local_file_name_gz,
                             local_file_name, bucket_name)
{
    setwd(local_base_folder_path)

    object_name <- file.path(remote_base_folder_path, file_name_gz)
    aws.s3::save_object(object = object_name,
                        key = account_key,
                        secret = account_secret,
                        bucket = bucket_name,
                        file = local_file_name_gz)

    #borra el gz (remove = FALSE ??):
    gunzip(local_file_name_gz, local_file_name)
}

upload_to_s3 <- function() {

    configure_aws_account()

    local_file_path <- '/u01/home/pbeati/.local/share/R/cladeAcc/100_way/raw/md5sum.txt'
    file.exists(local_file_path)

    # https://blog.djnavarro.net/posts/2022-03-17_using-aws-s3-in-r/

    bucket_exists(
        bucket = "cladeacc",
        region = "us-east-1"
    )

    bucket_exists("cladeacc")

    # download from s3

    tmp_storage_relative_path <- 'aws_s3'

    tmp_base_dir <- tools::R_user_dir("cladeAcc", which = "data")
    dir.create(tmp_base_dir, recursive = TRUE, showWarnings = FALSE)

    if(nchar(tmp_storage_relative_path) == 0) {
        data_base_dir <- tmp_base_dir
    } else {
        data_base_dir <- file.path(tmp_base_dir, tmp_storage_relative_path)
        dir.create(data_base_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # setwd(data_base_dir)

    out_file_path <- file.path(data_base_dir, 'md5sum.txt')

    #funciona
    # get_object("s3://cladeacc/100_way/md5sum.txt")

    #funciona
    # esto es download!!
    save_object(
        object = "100_way/md5sum.txt",
        bucket = "s3://cladeacc",
        region = "us-east-1",
        file = out_file_path
    )

    # upload to s3
    #funciona
    #upload
    put_object(
        file = local_file_path,
        object = '100_way/md5sum.txt',
        bucket = "cladeacc"
    )

    ?save_object
    ?put_object
}

# https://cran.r-project.org/web/packages/aws.s3/aws.s3.pdf
# https://www.rdocumentation.org/packages/aws.s3/versions/0.3.21

# bucketlist()
#
#
#
#
# ?aws.s3
#
# ?gunzip
#
# bucket_name <- 'cladeacc'
