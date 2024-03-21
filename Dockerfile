FROM janvaneck1994/msms_dssp

RUN pip3 install numpy
RUN pip3 install biopython scipy networkx pandas PyYAML

COPY ./MolPatch /MolPatch/

WORKDIR /MolPatch

ENTRYPOINT ["python3", "main.py"]
