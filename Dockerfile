FROM janvaneck1994/msms_dssp

RUN apt update && \
    apt install -y xvfb libxrender1 libxtst6 libxi6 libxkbcommon-x11-0 libxcomposite1

RUN pip3 install numpy biopython scipy networkx pandas PyYAML xvfbwrapper pyvirtualdisplay mayavi

COPY ./MolPatch /MolPatch/

WORKDIR /MolPatch

ENTRYPOINT ["python3", "main.py"]

