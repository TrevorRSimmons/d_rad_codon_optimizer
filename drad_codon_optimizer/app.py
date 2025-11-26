# app.py

import streamlit as st
from optimizer import optimize_sequence

st.set_page_config(page_title="The D. radiodurans Codon Optimizer", layout="centered")

st.markdown(
    """
    <div style="display:flex; align-items:center; justify-content:center; margin-bottom: 1rem;">
        <h1 style="margin:0; margin-right:15px;">
            The <em>Deinococcus radiodurans</em> Codon Optimizer
        </h1>
        <img src="https://upload.wikimedia.org/wikipedia/commons/7/73/Deinococcus_radiodurans._esquem%C3%A0tic.jpg"
             style="
                height:4em;
                width:4em;
                border-radius:50%;
                object-fit:cover;
                border: 2px solid #888;
             ">
    </div>
    """,
    unsafe_allow_html=True
)

st.markdown(
    """
    <p style='font-size:1.05em;'>
        Developed by
        <a href='https://scholar.google.com/citations?user=aSl7DJEAAAAJ&hl=en' target='_blank'>Laci Moline</a>,
        <a href='https://scholar.google.com/citations?user=pu-5jVYAAAAJ&hl=en' target='_blank'>Trevor Simmons Ph.D.</a>,
        and
        <a href='https://sites.utexas.edu/contreraslab/' target='_blank'>The Contreras Lab at UT-Austin</a>
        from the 
        <a href='https://pubs.acs.org/doi/10.1021/acssynbio.5c00409' target='_blank'>Standardized Genetic Toolkit for the Radiation-Resistant Extremophile Deinococcus radiodurans</a>,
    </p>
    """,
    unsafe_allow_html=True
)
st.markdown("Codon usage data sourced from the ""[Lowe Lab GtRNAdb dataset (Deinococcus radiodurans)](https://link.springer.com/protocol/10.1007/978-1-4939-9173-0_1)")

st.markdown(
    """
    <div style="margin-top:20px; margin-bottom:25px; font-size:1.05em;">
        Paste a nucleotide or amino acid sequence to obtain a <b>codon-optimized DNA</b>
        sequence for expression in <em>Deinococcus radiodurans</em>.
    </div>
    """,
    unsafe_allow_html=True
)

input_type_label = st.radio(
    "Input type:",
    ["Nucleotides", "Amino Acids"],
    horizontal=True,
)

seq = st.text_area(
    "Input sequence",
    height=200,
    placeholder="Paste your sequence here...",
)

wrap_output = st.checkbox("Wrap output at 60 nt/line", value=True)

if st.button("Optimize"):
    if not seq.strip():
        st.error("Please provide a sequence.")
    else:
        raw = seq.strip().replace(" ", "").replace("\n", "")

        if input_type_label == "Nucleotides":
            input_type = "dna"
            raw = raw.upper().replace("U", "T")
        else:
            input_type = "aa"
            raw = raw.upper()

        try:
            optimized = optimize_sequence(raw, input_type)

            if wrap_output:
                width = 60
                wrapped = "\n".join(
                    optimized[i:i+width] for i in range(0, len(optimized), width)
                )
            else:
                wrapped = optimized

            st.subheader("Optimized DNA sequence:")
            st.code(wrapped)

            st.info(f"Length: {len(optimized)} nucleotides")

            st.markdown("<p style='text-align:center; font-style:italic;'>Putting the rad in D. rad. since 2011</p>", unsafe_allow_html=True)
            
            st.download_button(
                "Download FASTA",
                data=f">optimized_Dradiodurans\n{wrapped}\n",
                file_name="optimized.fasta",
                mime="text/plain",
            )

        except ValueError as e:
            st.error(str(e))
        except Exception as e:
            st.error(f"Unexpected error: {e}")
