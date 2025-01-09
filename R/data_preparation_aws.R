library(aws.s3)
library(R.utils)

# Sys.setenv(
#     "AWS_DEFAULT_REGION" = "us-east-1"
# )

# "AWS_ACCESS_KEY_ID" = "mykey",
# "AWS_SECRET_ACCESS_KEY" = "mysecretkey",


# bucketlist(key = "s3_upload", secret = "", add_region = "us-east-1")




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



}

# https://cran.r-project.org/web/packages/aws.s3/aws.s3.pdf


# bucketlist()
#
# bucket_exists("cladeacc")
#
#
# ?aws.s3
#
# ?gunzip
#
# bucket_name <- 'cladeacc'
