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
        color: #2d3748;
    }
    
    .info-box {
        background: linear-gradient(135deg, #e0e7ff 0%, #f3e8ff 100%);
        border-left: 5px solid #667eea;
        border-radius: 10px;
        padding: 20px;
        margin: 15px 0;
        box-shadow: 0 5px 15px rgba(0,0,0,0.1);
    }
    
    .info-box p, .info-box div, .info-box span, .info-box li {
        color: #1a202c !important;
    }
    
    .info-box strong {
        color: #4c51bf !important;
        font-weight: 700 !important;
    }
    
    .metric-card {
        background: linear-gradient(135deg, #ffffff 0%, #f7fafc 100%);
        border-radius: 15px;
        padding: 20px;
        text-align: center;
        box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        transition: transform 0.3s ease;
        border: 2px solid #e2e8f0;
    }
    
    .metric-card:hover {
        transform: translateY(-5px);
    }
    
    .metric-card .metric-label {
        color: #718096;
        font-size: 0.9rem;
        font-weight: 600;
        margin-bottom: 8px;
    }
    
    .metric-card .metric-value {
        color: #2d3748;
        font-size: 1.5rem;
        font-weight: 700;
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
        background: linear-gradient(135deg, #fed7e2 0%, #feebc8 100%);
        color: #742a2a;
        padding: 15px;
        border-radius: 10px;
        margin: 10px 0;
        font-weight: 500;
        border-left: 5px solid #f56565;
    }
    
    .success-box {
        background: linear-gradient(135deg, #c6f6d5 0%, #9ae6b4 100%);
        color: #22543d;
        padding: 15px;
        border-radius: 10px;
        margin: 10px 0;
        font-weight: 500;
        border-left: 5px solid #48bb78;
    }
    
    .validation-highlight {
        background: #fef5e7;
        border: 2px solid #f39c12;
        border-radius: 8px;
        padding: 12px;
        margin: 10px 0;
        color: #7d6608 !important;
        font-weight: 600;
    }
    
    .stMarkdown p, .stMarkdown div, .stMarkdown span {
        color: #1a202c !important;
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
            coord_type = 'coding' if ':c.' in input_text else 'genomic'
            return {
                'accession': match.group(1),
                'position': int(match.group(2)),
                'ref': match.group(3).upper(),
                'alt': match.group(4).upper(),
                'type': coord_type
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

def fetch_genbank_record(accession):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        handle.close()
        return seq_record
    except Exception as e:
        st.error(f"Error fetching GenBank record: {str(e)}")
        return None

def extract_cds_sequence(genbank_record):
    for feature in genbank_record.features:
        if feature.type == "CDS":
            try:
                cds_seq = feature.extract(genbank_record.seq)
                return str(cds_seq), feature.location
            except:
                pass
    return None, None

def get_sequence_for_mutation(genbank_record, mutation_type):
    if mutation_type == 'coding':
        cds_seq, cds_location = extract_cds_sequence(genbank_record)
        if cds_seq:
            return cds_seq, 'CDS', cds_location
        else:
            st.warning("CDS not found in record, using full sequence")
            return str(genbank_record.seq), 'Full', None
    else:
        return str(genbank_record.seq), 'Genomic', None

def highlight_mutation(sequence, position, ref, alt):
    if position < 1 or position > len(sequence):
        return sequence, False, None
    
    idx = position - 1
    actual_base = sequence[idx].upper()
    
    if actual_base != ref.upper():
        highlighted = sequence[:idx] + f'<span class="mutation-highlight" style="background: #ffa500;">{actual_base}</span>' + sequence[idx+1:]
        return highlighted, False, actual_base
    
    highlighted = sequence[:idx] + f'<span class="mutation-highlight">{alt}</span>' + sequence[idx+1:]
    return highlighted, True, actual_base

def analyze_mutation(sequence, position, ref, alt):
    analysis = {}
    
    idx = position - 1
    
    if idx < 0 or idx >= len(sequence):
        analysis['validation'] = '❌ Position out of range'
        analysis['actual_base'] = 'N/A'
        return analysis
    
    actual_base = sequence[idx].upper()
    analysis['actual_base'] = actual_base
    
    if actual_base == ref.upper():
        analysis['validation'] = '✅ Reference allele matches'
        analysis['validation_status'] = 'match'
    else:
        analysis['validation'] = f'⚠️ Mismatch! Expected {ref}, found {actual_base}'
        analysis['validation_status'] = 'mismatch'
    
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
            analysis['effect'] = 'Synonymous'
        elif mutated_aa == '*':
            analysis['effect'] = 'Nonsense'
        elif original_aa == '*':
            analysis['effect'] = 'Readthrough'
        else:
            analysis['effect'] = 'Missense'
    
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
            placeholder="e.g., NM_022552.5:c.2688A>G or CDK9 c.208C>T",
            help="Supports multiple formats: NM_xxx:c.xxx (CDS), NC_xxx:g.xxx (genomic), GeneName c.xxx"
        )
    
    with col2:
        st.write("")
        st.write("")
        analyze_button = st.button("🔬 Analyze Mutation", width="stretch")
    
    st.markdown("<div class='info-box'>", unsafe_allow_html=True)
    st.markdown("""
    <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Supported Formats:</strong></p>
    <ul style='color: #1a202c !important;'>
    <li><code>NM_022552.5:c.2688A>G</code> - RefSeq with <strong>coding sequence</strong> position (CDS)</li>
    <li><code>NC_000009.11:g.130549830C>T</code> - Genomic position</li>
    <li><code>CDK9 c.208C>T</code> - Gene name with CDS mutation</li>
    <li><code>c.208C>T</code> - CDS position only (requires context)</li>
    </ul>
    <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Note:</strong> c. = coding sequence (CDS), g. = genomic coordinates</p>
    """, unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)
    
    if analyze_button and mutation_input:
        with st.spinner("🔍 Parsing mutation..."):
            time.sleep(0.3)
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
        
        with st.spinner(f"📥 Fetching GenBank record for {accession}..."):
            genbank_record = fetch_genbank_record(accession)
        
        if not genbank_record:
            st.error("Failed to fetch GenBank record from NCBI")
            return
        
        st.success(f"✅ Retrieved GenBank record: {len(genbank_record.seq)} bp")
        
        mutation_type = mutation_data.get('type', 'coding')
        
        sequence, seq_type, cds_location = get_sequence_for_mutation(genbank_record, mutation_type)
        
        if seq_type == 'CDS':
            st.info(f"📍 Working with CDS sequence: {len(sequence)} bp")
        elif seq_type == 'Genomic':
            st.info(f"📍 Working with genomic sequence: {len(sequence)} bp")
        else:
            st.info(f"📍 Working with full sequence: {len(sequence)} bp")
        
        position = mutation_data['position']
        ref = mutation_data['ref']
        alt = mutation_data['alt']
        
        highlighted_seq, is_valid, actual_base = highlight_mutation(sequence, position, ref, alt)
        
        analysis = analyze_mutation(sequence, position, ref, alt)
        
        if not is_valid:
            st.markdown(f"""<div class='validation-highlight'>
            ⚠️ <strong>VALIDATION WARNING:</strong> At position {position}, expected nucleotide <strong>{ref}</strong>, 
            but found <strong>{actual_base}</strong> in the sequence!
            </div>""", unsafe_allow_html=True)
        else:
            st.markdown(f"""<div class='success-box'>
            ✅ <strong>VALIDATION PASSED:</strong> Nucleotide at position {position} is <strong>{actual_base}</strong> (matches expected {ref})
            </div>""", unsafe_allow_html=True)
        
        st.markdown("---")
        st.subheader("📊 Mutation Analysis")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown("""
                <div class='metric-card'>
                    <div class='metric-label'>Position</div>
                    <div class='metric-value'>""" + str(position) + """</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown("""
                <div class='metric-card'>
                    <div class='metric-label'>Change</div>
                    <div class='metric-value'>""" + f"{ref}>{alt}" + """</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col3:
            st.markdown("""
                <div class='metric-card'>
                    <div class='metric-label'>Type</div>
                    <div class='metric-value'>""" + analysis.get('mutation_type', 'N/A') + """</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col4:
            st.markdown("""
                <div class='metric-card'>
                    <div class='metric-label'>Effect</div>
                    <div class='metric-value'>""" + analysis.get('effect', 'N/A') + """</div>
                </div>
            """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 🧪 Molecular Details")
            st.markdown("<div class='info-box'>", unsafe_allow_html=True)
            st.markdown(f"""
            <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Validation:</strong> <span style='color: #1a202c !important;'>{analysis['validation']}</span></p>
            <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Actual Base at Position:</strong> <span style='color: #1a202c !important; font-weight: 700;'>{analysis['actual_base']}</span></p>
            <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Sequence Type:</strong> <span style='color: #1a202c !important;'>{seq_type}</span></p>
            """, unsafe_allow_html=True)
            
            if 'original_codon' in analysis:
                st.markdown(f"""
                <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Original Codon:</strong> <span style='color: #1a202c !important;'>{analysis['original_codon']}</span></p>
                <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Mutated Codon:</strong> <span style='color: #1a202c !important;'>{analysis['mutated_codon']}</span></p>
                <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Codon Position:</strong> <span style='color: #1a202c !important;'>{analysis['codon_position']}</span></p>
                """, unsafe_allow_html=True)
            
            if 'original_aa' in analysis:
                st.markdown(f"""
                <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Amino Acid Change:</strong> <span style='color: #1a202c !important;'>{analysis['original_aa']} → {analysis['mutated_aa']}</span></p>
                """, unsafe_allow_html=True)
            
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("### 🎯 Sequence Context")
            st.markdown("<div class='info-box'>", unsafe_allow_html=True)
            st.markdown(f"""
            <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Accession:</strong> <span style='color: #1a202c !important;'>{accession}</span></p>
            <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Sequence Length:</strong> <span style='color: #1a202c !important;'>{len(sequence)} bp</span></p>
            """, unsafe_allow_html=True)
            
            if 'context' in analysis:
                st.markdown(f"""
                <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>Local Context:</strong> <code style='color: #1a202c !important; background: #f7fafc; padding: 2px 6px; border-radius: 3px;'>...{analysis['context']}...</code></p>
                """, unsafe_allow_html=True)
            
            if cds_location:
                st.markdown(f"""
                <p style='color: #1a202c !important;'><strong style='color: #4c51bf !important;'>CDS Location:</strong> <span style='color: #1a202c !important;'>{cds_location}</span></p>
                """, unsafe_allow_html=True)
            
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
            
            if is_valid:
                highlighted_base = f'<span class="mutation-highlight">{alt}</span>'
            else:
                highlighted_base = f'<span class="mutation-highlight" style="background: #ffa500;">{actual_base}</span>'
            
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
                width="stretch"
            )
        
        st.session_state.analysis_history.append({
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'input': mutation_input,
            'accession': accession,
            'position': position,
            'change': f"{ref}>{alt}",
            'effect': analysis.get('effect', 'N/A'),
            'validation': 'PASS' if is_valid else 'FAIL'
        })
    
    if st.session_state.analysis_history:
        st.markdown("---")
        st.subheader("📜 Analysis History")
        df = pd.DataFrame(st.session_state.analysis_history)
        st.dataframe(df, width="stretch", hide_index=True)

if __name__ == "__main__":
    main()
