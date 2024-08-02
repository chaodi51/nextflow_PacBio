process ALVIS {
    publishDir params.outdir, mode:'copy'
    publishDir params.datalakedir, mode:'copy'

    input:
        path paf

    output:
        path "chimeras.txt", emit: chimeras_txt

    script:
    """
        java -jar /home/circleci/project/alvis/dist/Alvis.jar -chimeras -printChimeras -chimeraPositions \
        -type contigalignment -inputfmt paf -outputfmt svg \
        -in ${paf} -outdir . -out test
        touch chimeras.txt
    """
}
