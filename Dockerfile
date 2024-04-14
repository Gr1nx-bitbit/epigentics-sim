FROM baseImage

WORKDIR /hypotheical/genomes

COPY . .

CMD [ "make", "run" ]