To load the galfast catalog of stars and the Opsim pointings history into a local postgresql database:

createdb calsim
gunzip -c selfcal.sql.gz | psql calsim

You can then change your environement variables:
LSST*_HOSTNAME to "localhost" 
LSST*_DATABASE to "calsim"
LSST*_USERNAME and LSST*_PASSWORD to "".
