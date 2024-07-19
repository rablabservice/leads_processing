import requests
import pandas as pd
from bs4 import BeautifulSoup

class FreeSurferColorLUT:
    def __init__(self):
        self.lut = self._fetch_lut()

    def _fetch_lut(self):
        url = "https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT"
        response = requests.get(url)
        soup = BeautifulSoup(response.content, "html.parser")
        pre_tag = soup.find("pre")
        lut_data = pre_tag.text
        lines = lut_data.split("\n")
        lut = []
        for line in lines:
            if not line.startswith("#") and line.strip() != "":
                parts = line.split()
                label_index = int(parts[0])
                label_name = " ".join(parts[1:-4])
                lut.append({"Label": label_index, "ROI": label_name})
        return lut

    def get_label_index_name(self, label_name):
        label_name = label_name.lower()
        label_index_name = []
        for label in self.lut:
            if label_name in label["ROI"].lower():
                label_index_name.append(label)
        if not label_index_name:
            raise ValueError(f"The ROI '{label_name}' does not exist in the FreeSurferColorLUT.\nFor more information, visit https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT")
        
        df = pd.DataFrame(label_index_name)
        return print(df.to_string(index=False))

    # Get just the label name from the label index

    def get_label_name(self, label_index):
        label_name = []
        for label in self.lut:
            if label_index == label["Label"]:
                label_name.append(label)
        if not label_name:
            raise ValueError(f"The label index '{label_index}' does not exist in the FreeSurferColorLUT.\nFor more information, visit https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT")

        return label_name[0]['ROI']
