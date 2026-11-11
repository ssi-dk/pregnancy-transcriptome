COMMANDS MADE TO RUN DECONVOLUTION ANALYSIS

docker pull cibersortx/fractions
docker images # just to check that I have succesfully obtained the image

COMMAND TO RUN CIBERSORTx:
docker run -v absolute/path/to/input/dir:/src/data -v absolute/path/to/output/dir:/src/outdir cibersortx/fractions --username email_address_registered_on_CIBERSORTx_website --token token_obtained_from_CIBERSORTx_website --mixture mixture_file.txt --sigmatrix LM22.txt --perm 100 --rmbatchBmode TRUE

ADD CIBERSORT VERSION INFO
