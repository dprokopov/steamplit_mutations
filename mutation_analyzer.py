import streamlit as st
import re
from Bio import Entrez, SeqIO
from io import StringIO
import time
import pandas as pd
from datetime import datetime

st.set_page_config(page_title="Mutation Analyzer", page_icon="🧬", layout="wide")

Entrez.email = "mutation.analyzer@example.com"

def apply_custom_css():
    st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;600;700&display=swap');
    
    * {
        font-family: 'Inter', sans-serif;
    }
    
    .main {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        animation: gradientShift 15s ease infinite;
    }
    
    @keyframes gradientShift {
        0%, 100% { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); }
        50% { background: linear-gradient(135deg, #764ba2 0%, #667eea 100%); }
    }
    
    .stApp {
        background: transparent;
    }
    
    div[data-testid="stVerticalBlock"] > div:has(div.element-container) {
        background: rgba(255, 255, 255, 0.95);
        backdrop-filter: blur(10px);
        border-radius: 20px;
        padding: 30px;
        box-shadow: 0 20px 60px rgba(0,0,0,0.3);
        margin: 20px 0;
        animation: fadeInUp 0.6s ease;
    }
    
    @keyframes fadeInUp {
        from {
            opacity: 0;
            transform: translateY(30px);
        }
        to {
            opacity: 1;
            transform: translateY(0);
        }
    }
    
    h1 {
        background: linear-gradient(120deg, #667eea, #764ba2);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-weight: 700;
        font-size: 3rem !important;
        text-align: center;
        margin-bottom: 10px;
        animation: titlePulse 3s ease-in-out infinite;
    }
    
    @keyframes titlePulse {
        0%, 100% { opacity: 1; }
        50% { opacity: 0.8; }
    }
    
    .subtitle {
        text-align: center;
        color: #ffffff;
        font-size: 1.2rem;
        margin-bottom: 30px;
        text-shadow: 2px 2px 4px rgba(0,0,0,0.2);
    }
    
    .stTextInput input {
        border-radius: 15px;
        border: 2px solid #667eea;
        padding: 15px;
        font-size: 1rem;
        transition: all 0.3s ease;
    }
    
    .stTextInput input:focus {
        border-color: #764ba2;
        box-shadow: 0 0 20px rgba(118, 75, 162, 0.3);
        transform: scale(1.02);
    }
    
    .stButton button {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 15px;
        padding: 15px 40px;
        font-size: 1.1rem;
        font-weight: 600;
        transition: all 0.3s ease;
        box-shadow: 0 10px 30px rgba(102, 126, 234, 0.4);
    }
    
    .stButton button:hover {
        transform: translateY(-3px);
        box-shadow: 0 15px 40px rgba(102, 126, 234, 0.6);
    }
    
    .mutation-highlight {
        background: linear-gradient(120deg, #f093fb 0%, #f5576c 100%);
        color: white;
        padding: 2px 6px;
        border-radius: 5px;
        font-weight: bold;
        animation: glow 2s ease-in-out infinite;
    }
    
    @keyframes glow {
        0%, 100% { box-shadow: 0 0 10px rgba(245, 87, 108, 0.5); }
        50% { box-shadow: 0 0 20px rgba(245, 87, 108, 0.8); }
    }
    
    .sequence-container {
        background: #f8f9fa;
        border-radius: 15px;
        padding: 20px;
        font-family: 'Courier New', monospace;
        font-size: 14px;
        line-height: 1.8;
        word-wrap: break-word;
        box-shadow: inset 0 2px 10px rgba(0,0,0,0.1);
        max-height: 400px;
        overflow-y: auto;
    }
    
    .info-box {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
        border-left: 5px solid #667eea;
        border-radius: 10px;
        padding: 20px;
        margin: 15px 0;
        box-shadow: 0 5px 15px rgba(0,0,0,0.1);
    }
    
    .metric-card {
        background: white;
        border-radius: 15px;
        padding: 20px;
        text-align: center;
        box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        transition: transform 0.3s ease;
    }
    
    .metric-card:hover {
        transform: translateY(-5px);
    }
    
    .stDownloadButton button {
        background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%);
        color: white;
        border: none;
        border-radius: 15px;
        padding: 12px 30px;
        font-weight: 600;
    }
    
    .warning-box {
        background: linear-gradient(135deg, #fa709a 0%, #fee140 100%);
        color: white;
        padding: 15px;
        border-radius: 10px;
        margin: 10px 0;
        font-weight: 500;
    }
    
    .success-box {
        background: linear-gradient(135deg, #13f1fc 0%, #0470dc 100%);
        color: white;
        padding: 15px;
        border-radius: 10px;
        margin: 10px 0;
        font-weight: 500;
    }
    </style>
    """, unsafe_allow_html=True)

def parse_mutation_input(input_text):
    input_text = input_text.strip().replace(" ", "")
    
    patterns = [
        r'(NM_\d+\.\d+):c\.(\d+)([ACGT])>([ACGT])',
        r'(NC_\d+\.\d+):g\.(\d+)([ACGT])>([ACGT])',
        r'(NM_\d+):c\.(\d+)([ACGT])>([ACGT])',
        r'(NC_\d+):g\.(\d+)([ACGT])>([ACGT])',
    ]
    
    for pattern in patterns:
        match = re.search(pattern, input_text, re.IGNORECASE)
        if match:
            return {
                'accession': match.group(1),
                'position': int(match.group(2)),
                'ref': match.group(3).upper(),
                'alt': match.group(4).upper(),
                'type': 'coding' if 'c.' in input_text else 'genomic'
            }
    
    gene_pattern = r'([A-Z0-9]+).*c\.(\d+)([ACGT])>([ACGT])'
    match = re.search(gene_pattern, input_text, re.IGNORECASE)
    if match:
        return {
            'gene': match.group(1).upper(),
            'position': int(match.group(2)),
            'ref': match.group(3).upper(),
            'alt': match.group(4).upper(),
            'type': 'gene_based'
        }
    
    simple_pattern = r'c\.(\d+)([ACGT])>([ACGT])'
    match = re.search(simple_pattern, input_text, re.IGNORECASE)
    if match:
        return {
            'position': int(match.group(1)),
            'ref': match.group(2).upper(),
            'alt': match.group(3).upper(),
            'type': 'position_only'
        }
    
    return None

def search_gene_accession(gene_name):
    try:
        search_term = f"{gene_name}[Gene] AND human[Organism] AND refseq[Filter]"
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        
        if record['IdList']:
            for seq_id in record['IdList']:
                handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
                seq_record = SeqIO.read(handle, "genbank")
                handle.close()
                
                if seq_record.id.startswith('NM_'):
                    return seq_record.id
            
            handle = Entrez.efetch(db="nucleotide", id=record['IdList'][0], rettype="gb", retmode="text")
            seq_record = SeqIO.read(handle, "genbank")
            handle.close()
            return seq_record.id
    except:
        pass
    return None

def fetch_sequence(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(seq_record.seq)
    except Exception as e:
        st.error(f"Error fetching sequence: {str(e)}")
        return None

def highlight_mutation(sequence, position, ref, alt):
    if position < 1 or position > len(sequence):
        return sequence, False
    
    idx = position - 1
    if sequence[idx].upper() != ref.upper():
        return sequence, False
    
    highlighted = sequence[:idx] + f'<span class="mutation-highlight">{alt}</span>' + sequence[idx+1:]
    return highlighted, True

def analyze_mutation(sequence, position, ref, alt):
    analysis = {}
    
    idx = position - 1
    
    if sequence[idx].upper() == ref.upper():
        analysis['validation'] = '✅ Reference allele matches'
    else:
        analysis['validation'] = f'⚠️ Reference mismatch: expected {ref}, found {sequence[idx]}'
    
    if ref != alt:
        if (ref in 'AG' and alt in 'AG') or (ref in 'CT' and alt in 'CT'):
            analysis['mutation_type'] = 'Transition'
        else:
            analysis['mutation_type'] = 'Transversion'
    
    codon_start = (idx // 3) * 3
    codon_end = codon_start + 3
    
    if codon_end <= len(sequence):
        original_codon = sequence[codon_start:codon_end]
        position_in_codon = idx % 3
        mutated_codon = original_codon[:position_in_codon] + alt + original_codon[position_in_codon+1:]
        
        analysis['original_codon'] = original_codon.upper()
        analysis['mutated_codon'] = mutated_codon.upper()
        analysis['codon_position'] = position_in_codon + 1
        
        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        original_aa = codon_table.get(original_codon.upper(), '?')
        mutated_aa = codon_table.get(mutated_codon.upper(), '?')
        
        analysis['original_aa'] = original_aa
        analysis['mutated_aa'] = mutated_aa
        
        if original_aa == mutated_aa:
            analysis['effect'] = 'Synonymous (Silent)'
        elif mutated_aa == '*':
            analysis['effect'] = 'Nonsense (Stop gain)'
        elif original_aa == '*':
            analysis['effect'] = 'Readthrough (Stop loss)'
        else:
            analysis['effect'] = 'Missense (Amino acid change)'
    
    context_start = max(0, idx - 10)
    context_end = min(len(sequence), idx + 11)
    analysis['context'] = sequence[context_start:context_end]
    
    return analysis

def create_download_content(sequence, position, ref, alt, accession):
    idx = position - 1
    marked_sequence = sequence[:idx] + f'*{alt}*' + sequence[idx+1:]
    
    content = f">Mutated_{accession}_pos{position}_{ref}>{alt}\n"
    
    for i in range(0, len(marked_sequence), 60):
        content += marked_sequence[i:i+60] + "\n"
    
    return content

def main():
    apply_custom_css()
    
    st.markdown("<h1>🧬 Mutation Analyzer</h1>", unsafe_allow_html=True)
    st.markdown("<p class='subtitle'>Advanced DNA Mutation Analysis & Visualization Platform</p>", unsafe_allow_html=True)
    
    if 'analysis_history' not in st.session_state:
        st.session_state.analysis_history = []
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        mutation_input = st.text_input(
            "Enter Mutation Identifier",
            placeholder="e.g., NM_001261.4:c.208C>T or CDK9 c.208C>T",
            help="Supports multiple formats: NM_xxx:c.xxx, NC_xxx:g.xxx, GeneName c.xxx"
        )
    
    with col2:
        st.write("")
        st.write("")
        analyze_button = st.button("🔬 Analyze Mutation", use_container_width=True)
    
    st.markdown("<div class='info-box'>", unsafe_allow_html=True)
    st.markdown("""
    **Supported Formats:**
    - `NM_001261.4:c.208C>T` - RefSeq with coding position
    - `NC_000009.11:g.130549830C>T` - Genomic position
    - `CDK9 c.208C>T` - Gene name with mutation
    - `c.208C>T` - Position only (requires context)
    """)
    st.markdown("</div>", unsafe_allow_html=True)
    
    if analyze_button and mutation_input:
        with st.spinner("🔍 Parsing mutation..."):
            time.sleep(0.5)
            mutation_data = parse_mutation_input(mutation_input)
        
        if not mutation_data:
            st.markdown("<div class='warning-box'>⚠️ Unable to parse mutation format. Please check your input.</div>", unsafe_allow_html=True)
            return
        
        st.markdown("<div class='success-box'>✅ Mutation parsed successfully!</div>", unsafe_allow_html=True)
        
        accession = None
        if 'accession' in mutation_data:
            accession = mutation_data['accession']
        elif 'gene' in mutation_data:
            with st.spinner(f"🔎 Searching for {mutation_data['gene']} gene..."):
                accession = search_gene_accession(mutation_data['gene'])
            
            if accession:
                st.success(f"Found accession: {accession}")
            else:
                st.error(f"Could not find RefSeq accession for gene {mutation_data['gene']}")
                return
        else:
            st.warning("No accession or gene provided. Please specify a complete mutation identifier.")
            return
        
        with st.spinner(f"📥 Fetching sequence for {accession}..."):
            sequence = fetch_sequence(accession)
        
        if not sequence:
            st.error("Failed to fetch sequence from NCBI")
            return
        
        st.success(f"✅ Retrieved sequence: {len(sequence)} bp")
        
        position = mutation_data['position']
        ref = mutation_data['ref']
        alt = mutation_data['alt']
        
        highlighted_seq, is_valid = highlight_mutation(sequence, position, ref, alt)
        
        if not is_valid:
            st.warning(f"⚠️ Position {position} may be out of range or reference allele doesn't match")
        
        st.markdown("---")
        st.subheader("📊 Mutation Analysis")
        
        analysis = analyze_mutation(sequence, position, ref, alt)
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown("<div class='metric-card'>", unsafe_allow_html=True)
            st.metric("Position", f"{position}")
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("<div class='metric-card'>", unsafe_allow_html=True)
            st.metric("Change", f"{ref}>{alt}")
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col3:
            st.markdown("<div class='metric-card'>", unsafe_allow_html=True)
            st.metric("Type", analysis.get('mutation_type', 'N/A'))
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col4:
            st.markdown("<div class='metric-card'>", unsafe_allow_html=True)
            st.metric("Effect", analysis.get('effect', 'N/A'))
            st.markdown("</div>", unsafe_allow_html=True)
        
        st.markdown("---")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 🧪 Molecular Details")
            st.markdown("<div class='info-box'>", unsafe_allow_html=True)
            st.write(f"**Validation:** {analysis['validation']}")
            if 'original_codon' in analysis:
                st.write(f"**Original Codon:** {analysis['original_codon']}")
                st.write(f"**Mutated Codon:** {analysis['mutated_codon']}")
                st.write(f"**Codon Position:** {analysis['codon_position']}")
            if 'original_aa' in analysis:
                st.write(f"**Amino Acid Change:** {analysis['original_aa']} → {analysis['mutated_aa']}")
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("### 🎯 Sequence Context")
            st.markdown("<div class='info-box'>", unsafe_allow_html=True)
            st.write(f"**Accession:** {accession}")
            st.write(f"**Sequence Length:** {len(sequence)} bp")
            if 'context' in analysis:
                st.write(f"**Local Context:** ...{analysis['context']}...")
            st.markdown("</div>", unsafe_allow_html=True)
        
        st.markdown("---")
        st.subheader("🧬 Sequence Visualization")
        
        show_full = st.checkbox("Show full sequence", value=False)
        
        if show_full:
            display_seq = highlighted_seq
        else:
            context_range = 200
            idx = position - 1
            start = max(0, idx - context_range)
            end = min(len(sequence), idx + context_range)
            
            pre_context = sequence[start:idx]
            post_context = sequence[idx+1:end]
            highlighted_base = f'<span class="mutation-highlight">{alt}</span>'
            display_seq = pre_context + highlighted_base + post_context
        
        st.markdown(f"<div class='sequence-container'>{display_seq}</div>", unsafe_allow_html=True)
        
        st.markdown("---")
        
        download_content = create_download_content(sequence, position, ref, alt, accession)
        
        col1, col2, col3 = st.columns([1, 1, 1])
        
        with col2:
            st.download_button(
                label="⬇️ Download Mutated Sequence (FASTA)",
                data=download_content,
                file_name=f"mutated_{accession}_{position}_{ref}_{alt}.fasta",
                mime="text/plain",
                use_container_width=True
            )
        
        st.session_state.analysis_history.append({
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'input': mutation_input,
            'accession': accession,
            'position': position,
            'change': f"{ref}>{alt}",
            'effect': analysis.get('effect', 'N/A')
        })
    
    if st.session_state.analysis_history:
        st.markdown("---")
        st.subheader("📜 Analysis History")
        df = pd.DataFrame(st.session_state.analysis_history)
        st.dataframe(df, use_container_width=True, hide_index=True)

if __name__ == "__main__":
    main()
