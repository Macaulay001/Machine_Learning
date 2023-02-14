import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle

st.set_page_config(
   page_title="Acetylcholine",
   page_icon="🧊",
   layout="wide",
#    initial_sidebar_state="expanded",
)
# Molecular descriptor calculator
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def acetylcholinesterase_build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

def glucagon_build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('glucagon_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

def gcprs_build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('gcprs_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

# Logo image
image = Image.open('logo.png')

st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Compounds Bioactivity Prediction App for Drug Target in Metabolic Syndrome

This app allows you to predict the bioactivity towards inhibting the `Acetylcholinesterase` enzyme. `Acetylcholinesterase` is a drug target for Alzheimer's disease.

**Credits**
- App built in `Python` + `Streamlit` by Macaulay Oladimeji (aka Textbook)
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
- Data Professor
---
""")

# Sidebar
st.header('Upload your CSV data')
uploaded_file = st.file_uploader("Upload your input file", type=['txt'])
st.markdown("""
[Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
""")

receptor_mode = st.selectbox('Choose a Target Receptor',['Acetylcholinesterase', 'Glucagon Receptor', 'GCPRs'])







if st.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)

    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    desc = pd.read_csv('descriptors_output.csv')
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    if receptor_mode=='Acetylcholinesterase':
        st.header('**Subset of descriptors from previously built Acetylcholinesterase models**')
        Xlist = list(pd.read_csv('acetylcholine_descriptor_list.csv').columns)
        desc_subset = desc[Xlist]
        st.write(desc_subset)
        st.write(desc_subset.shape)
        # Apply trained model to make prediction on query compounds
        acetylcholinesterase_build_model(desc_subset)
    elif receptor_mode=='Glucagon Receptor':
        st.header('**Subset of descriptors from previously built Glucagon Receptor models**')
        Xlist = list(pd.read_csv('glucagon_descriptor_list.csv').columns)
        desc_subset = desc[Xlist]
        st.write(desc_subset)
        st.write(desc_subset.shape)
        # Apply trained model to make prediction on query compounds
        glucagon_build_model(desc_subset)
    elif receptor_mode=='GCPRs':
        st.header('**Subset of descriptors from previously built GCPRs models**')
        Xlist = list(pd.read_csv('gcprs_descriptor_list.csv').columns)
        desc_subset = desc[Xlist]
        st.write(desc_subset)
        st.write(desc_subset.shape)
        # Apply trained model to make prediction on query compounds
        gcprs_build_model(desc_subset)
    else:
        st.info('Reload Page and choose a target Receptor!')

    
else:
    st.info('Upload input data to start!')
