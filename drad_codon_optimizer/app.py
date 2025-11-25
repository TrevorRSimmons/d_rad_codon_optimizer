# app.py

import streamlit as st
from optimizer import optimize_sequence

st.set_page_config(page_title="The D. radiodurans Codon Optimizer", layout="centered")

st.title("The Deinococcus radiodurans Codon Optimizer")
st.caption("Originally developed by Laci Moline, Trevor Simmons, and The Contreras Lab at UT-Austin")

st.markdown("""
Paste a nucleotide or amino acid sequence to obtain a **codon-optimized DNA**
sequence for expression in *Deinococcus radiodurans*.
""")

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

            st.markdown("<p style='text-align:center; font-style:italic;'>Putting the rad in D. rad. since 2011!</p>", unsafe_allow_html=True)
            
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
