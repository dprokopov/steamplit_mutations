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
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    
    * {
        font-family: 'Inter', sans-serif;
    }
    
    .main {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    }
    
    .stApp {
        background: transparent;
    }
    
    div[data-testid="stVerticalBlock"] > div:has(div.element-container) {
        background: rgba(255, 255, 255, 0.98);
        backdrop-filter: blur(10px);
        border-radius: 16px;
        padding: 32px;
        box-shadow: 0 8px 32px rgba(0,0,0,0.1);
        margin: 20px 0;
    }
    
    h1 {
        background: white;
        color: #667eea;
        font-weight: 700;
        font-size: 2.5rem !important;
        text-align: center;
        margin-bottom: 16px;
        padding: 24px;
        border-radius: 12px;
        box-shadow: 0 4px 20px rgba(0,0,0,0.1);
    }
    
    h2, h3 {
        color: #2d3748 !important;
        font-weight: 600 !important;
        margin-top: 24px !important;
        margin-bottom: 16px !important;
    }
    
    .subtitle {
        text-align: center;
        color: #2d3748 !important;
        font-size: 1rem;
        margin-bottom: 24px;
        background: rgba(255, 255, 255, 0.95);
        padding: 12px 24px;
        border-radius: 8px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        font-weight: 500;
    }
    
    .stTextInput input {
        border-radius: 8px;
        border: 2px solid #e2e8f0;
        padding: 12px 16px;
        font-size: 0.95rem;
        transition: all 0.2s ease;
        color: #2d3748;
    }
    
    .stTextInput input:focus {
        border-color: #667eea;
        box-shadow: 0 0 0 3px rgba(102, 126, 234, 0.1);
    }
    
    .stButton button {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 12px 32px;
        font-size: 1rem;
        font-weight: 600;
        transition: all 0.2s ease;
        box-shadow: 0 4px 12px rgba(102, 126, 234, 0.3);
    }
    
    .stButton button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 16px rgba(102, 126, 234, 0.4);
    }
    
    .mutation-highlight {
        background: #f56565;
        color: white;
        padding: 2px 6px;
        border-radius: 4px;
        font-weight: 600;
    }
    
    .sequence-container {
        background: #f7fafc;
        border-radius: 8px;
        padding: 20px;
        font-family: 'Courier New', monospace;
        font-size: 14px;
        line-height: 1.8;
        word-wrap: break-word;
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.06);
        max-height: 400px;
        overflow-y: auto;
        color: #2d3748;
        border: 1px solid #e2e8f0;
    }
    
    .info-box {
        background: white;
        border-left: 4px solid #667eea;
        border-radius: 8px;
        padding: 20px;
        margin: 16px 0;
        box-shadow: 0 2px 8px rgba(0,0,0,0.05);
    }
    
    .info-box p, .info-box div, .info-box span, .info-box li, .info-box ul, .info-box code {
        color: #2d3748 !important;
        font-size: 0.95rem;
        line-height: 1.6;
    }
    
    .info-box strong {
        color: #667eea !important;
        font-weight: 600 !important;
    }
    
    .info-box code {
        background: #f7fafc !important;
        padding: 2px 6px !important;
        border-radius: 4px !important;
        border: 1px solid #e2e8f0 !important;
        font-family: 'Courier New', monospace !important;
        font-size: 0.9rem !important;
    }
    
    .info-box ul {
        margin: 8px 0;
        padding-left: 24px;
    }
    
    .info-box li {
        margin: 4px 0;
    }
    
    .metric-card {
        background: white;
        border-radius: 8px;
        padding: 20px;
        text-align: center;
        box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        transition: transform 0.2s ease;
        border: 1px solid #e2e8f0;
    }
    
    .metric-card:hover {
        transform: translateY(-2px);
        border-color: #667eea;
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
    }
    
    .metric-card .metric-label {
        color: #718096;
        font-size: 0.85rem;
        font-weight: 600;
        margin-bottom: 8px;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    
    .metric-card .metric-value {
        color: #667eea;
        font-size: 1.5rem;
        font-weight: 700;
    }
    
    .stDownloadButton button {
        background: linear-gradient(135deg, #48bb78 0%, #38a169 100%);
        color: white;
        border: none;
        border-radius: 8px;
        padding: 10px 24px;
        font-weight: 600;
        box-shadow: 0 2px 8px rgba(72, 187, 120, 0.3);
    }
    
    .stDownloadButton button:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(72, 187, 120, 0.4);
    }
    
    .warning-box {
        background: #fff5f5;
        color: #c53030;
        padding: 16px;
        border-radius: 8px;
        margin: 16px 0;
        font-weight: 500;
        border: 2px solid #fc8181;
        font-size: 0.95rem;
    }
    
    .success-box {
        background: #f0fff4;
        color: #22543d;
        padding: 16px;
        border-radius: 8px;
        margin: 16px 0;
        font-weight: 500;
        border: 2px solid #48bb78;
        font-size: 0.95rem;
    }
    
    .validation-highlight {
        background: #fffaf0;
        border: 2px solid #ed8936;
        border-radius: 8px;
        padding: 16px;
        margin: 16px 0;
        color: #7c2d12 !important;
        font-weight: 500;
        font-size: 0.95rem;
    }
    
    .stMarkdown {
        color: #2d3748 !important;
    }
    
    .stMarkdown p, .stMarkdown div, .stMarkdown span {
        color: #2d3748 !important;
    }
    
    .stSpinner > div {
        border-color: #667eea !important;
    }
    
    hr {
        margin: 32px 0;
        border: none;
        border-top: 2px solid #e2e8f0;
    }
    
    .stAlert {
        border-radius: 8px;
        padding: 12px 16px;
    }
    </style>
    """, unsafe_allow_html=True)

def parse_mutation_input(input_text):
    input_text = input_text.strip()
    
    simple_sub_pattern = r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):c\.(\d+)([ACGT])>([ACGT])'
    match = re.search(simple_sub_pattern, input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'position': int(match.group(3)),
            'ref': match.group(4).upper(),
            'alt': match.group(5).upper(),
            'type': 'coding'
        }
    
    intronic_pattern = r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):c\.(\d+)([\+\-])(\d+)([ACGT])>([ACGT])'
    match = re.search(intronic_pattern, input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'position': int(match.group(3)),
            'intronic_offset': f"{match.group(4)}{match.group(5)}",
            'ref': match.group(6).upper(),
            'alt': match.group(7).upper(),
            'type': 'intronic'
        }
    
    deletion_pattern = r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):c\.(\d+)([\+\-])?(\d+)?del'
    match = re.search(deletion_pattern, input_text, re.IGNORECASE)
    if match:
        offset = f"{match.group(4)}{match.group(5)}" if match.group(4) else None
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'position': int(match.group(3)),
            'intronic_offset': offset,
            'type': 'deletion'
        }
    
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

def get_sequence_safely(seq_obj):
    try:
        return str(seq_obj)
    except Exception:
        return None

def get_sequence_for_mutation(genbank_record, mutation_type, position=None):
    if mutation_type == 'coding':
        cds_seq, cds_location = extract_cds_sequence(genbank_record)
        if cds_seq:
            return cds_seq, 'CDS', cds_location
        else:
            st.warning("CDS not found in record, using full sequence")
    
    seq_length = len(genbank_record.seq)
    
    if seq_length > 10000000:
        st.warning(f"⚠️ Large genomic sequence detected ({seq_length:,} bp). Extracting relevant region...")
        
        if position:
            context = 5000
            start = max(0, position - context)
            end = min(seq_length, position + context)
            
            try:
                sub_seq = genbank_record.seq[start:end]
                seq_str = get_sequence_safely(sub_seq)
                
                if seq_str:
                    adjusted_position = position - start
                    return seq_str, 'Genomic (Region)', None, adjusted_position, start
                else:
                    st.error("Unable to extract sequence from large genomic record")
                    return None, None, None, None, None
            except Exception as e:
                st.error(f"Error extracting region: {str(e)}")
                return None, None, None, None, None
        else:
            st.error("Position required to extract from large genomic sequence")
            return None, None, None, None, None
    
    seq_str = get_sequence_safely(genbank_record.seq)
    if seq_str:
        return seq_str, 'Genomic', None, None, None
    else:
        st.error("Unable to read sequence from record")
        return None, None, None, None, None

def highlight_mutation(sequence, position, ref, alt):
    if position < 1 or position > len(sequence):
        return sequence, False, None
    
    idx = position - 1
    actual_base = sequence[idx].upper()
    
    if actual_base != ref.upper():
        highlighted = sequence[:idx] + f'<span class="mutation-highlight" style="background: #ed8936;">{actual_base}</span>' + sequence[idx+1:]
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
    st.markdown("<p class='subtitle'>Advanced genomic mutation analysis tool</p>", unsafe_allow_html=True)
    
    if 'analysis_history' not in st.session_state:
        st.session_state.analysis_history = []
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        mutation_input = st.text_input(
            "Enter Mutation Identifier",
            placeholder="e.g., NM_002440.4(MSH4):c.23C>T or NM_022552.5:c.2688A>G",
            help="Supports RefSeq with CDS positions, intronic variants, and deletions"
        )
    
    with col2:
        st.write("")
        st.write("")
        analyze_button = st.button("🔬 Analyze", use_container_width=True)
    
    with st.expander("ℹ️ Supported Formats", expanded=False):
        st.markdown("""
        **Coding Sequence Variants:**
        - `NM_002440.4(MSH4):c.23C>T` - RefSeq with gene name
        - `NM_022552.5:c.2688A>G` - RefSeq with CDS position
        
        **Intronic Variants:**
        - `NM_002440.4(MSH4):c.1906+8G>T` - Intronic position
        
        **Deletions:**
        - `NM_002440.4(MSH4):c.2620-10del` - Deletion variant
        
        **Note:** For best results, use NM_ RefSeq accessions with CDS coordinates (c.)
        """)
    
    if analyze_button and mutation_input:
        with st.spinner("🔍 Parsing mutation..."):
            time.sleep(0.2)
            mutation_data = parse_mutation_input(mutation_input)
        
        if not mutation_data:
            st.markdown("<div class='warning-box'>⚠️ Unable to parse mutation format. Please check the supported formats above.</div>", unsafe_allow_html=True)
            return
        
        if mutation_data['type'] in ['deletion', 'intronic']:
            st.markdown(f"<div class='warning-box'>⚠️ {mutation_data['type'].capitalize()} variants are detected but not fully supported for analysis yet. Only simple substitutions are currently analyzed.</div>", unsafe_allow_html=True)
            return
        
        st.markdown("<div class='success-box'>✅ Mutation parsed successfully</div>", unsafe_allow_html=True)
        
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
        
        st.success(f"✅ Retrieved record: {len(genbank_record.seq):,} bp")
        
        mutation_type = mutation_data.get('type', 'coding')
        position = mutation_data['position']
        
        result = get_sequence_for_mutation(genbank_record, mutation_type, position)
        
        if result[0] is None:
            st.error("Failed to extract sequence")
            return
        
        sequence = result[0]
        seq_type = result[1]
        cds_location = result[2]
        adjusted_position = result[3] if len(result) > 3 and result[3] else position
        region_start = result[4] if len(result) > 4 else None
        
        if seq_type == 'CDS':
            st.info(f"📍 Working with CDS sequence: {len(sequence):,} bp")
        elif 'Region' in seq_type:
            st.info(f"📍 Extracted genomic region: {len(sequence):,} bp (Position adjusted within region)")
        else:
            st.info(f"📍 Working with {seq_type.lower()} sequence: {len(sequence):,} bp")
        
        ref = mutation_data['ref']
        alt = mutation_data['alt']
        
        highlighted_seq, is_valid, actual_base = highlight_mutation(sequence, adjusted_position, ref, alt)
        
        analysis = analyze_mutation(sequence, adjusted_position, ref, alt)
        
        if not is_valid:
            st.markdown(f"""<div class='validation-highlight'>
            ⚠️ <strong>VALIDATION WARNING:</strong> At position {position}, expected <strong>{ref}</strong>, 
            but found <strong>{actual_base}</strong> in the sequence
            </div>""", unsafe_allow_html=True)
        else:
            st.markdown(f"""<div class='success-box'>
            ✅ <strong>VALIDATED:</strong> Nucleotide at position {position} is <strong>{actual_base}</strong> (matches {ref})
            </div>""", unsafe_allow_html=True)
        
        st.markdown("---")
        st.subheader("📊 Mutation Analysis")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Position</div>
                    <div class='metric-value'>{position}</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col2:
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Change</div>
                    <div class='metric-value'>{ref}>{alt}</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col3:
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Type</div>
                    <div class='metric-value'>{analysis.get('mutation_type', 'N/A')}</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col4:
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Effect</div>
                    <div class='metric-value'>{analysis.get('effect', 'N/A')}</div>
                </div>
            """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 🧪 Molecular Details")
            st.markdown("<div class='info-box'>", unsafe_allow_html=True)
            st.markdown(f"<p><strong>Validation:</strong> {analysis['validation']}</p>", unsafe_allow_html=True)
            st.markdown(f"<p><strong>Actual Base:</strong> {analysis['actual_base']}</p>", unsafe_allow_html=True)
            st.markdown(f"<p><strong>Sequence Type:</strong> {seq_type}</p>", unsafe_allow_html=True)
            
            if 'original_codon' in analysis:
                st.markdown(f"<p><strong>Original Codon:</strong> {analysis['original_codon']}</p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Mutated Codon:</strong> {analysis['mutated_codon']}</p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Codon Position:</strong> {analysis['codon_position']}</p>", unsafe_allow_html=True)
            
            if 'original_aa' in analysis:
                st.markdown(f"<p><strong>Amino Acid Change:</strong> {analysis['original_aa']} → {analysis['mutated_aa']}</p>", unsafe_allow_html=True)
            
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("### 🎯 Sequence Context")
            st.markdown("<div class='info-box'>", unsafe_allow_html=True)
            st.markdown(f"<p><strong>Accession:</strong> {accession}</p>", unsafe_allow_html=True)
            st.markdown(f"<p><strong>Sequence Length:</strong> {len(sequence):,} bp</p>", unsafe_allow_html=True)
            
            if 'context' in analysis:
                st.markdown(f"<p><strong>Local Context:</strong> <code>...{analysis['context']}...</code></p>", unsafe_allow_html=True)
            
            if cds_location:
                st.markdown(f"<p><strong>CDS Location:</strong> {cds_location}</p>", unsafe_allow_html=True)
            
            if region_start:
                st.markdown(f"<p><strong>Genomic Start:</strong> {region_start:,}</p>", unsafe_allow_html=True)
            
            st.markdown("</div>", unsafe_allow_html=True)
        
        st.markdown("---")
        st.subheader("🧬 Sequence Visualization")
        
        show_full = st.checkbox("Show full sequence", value=False)
        
        if show_full:
            display_seq = highlighted_seq
        else:
            context_range = 150
            idx = adjusted_position - 1
            start = max(0, idx - context_range)
            end = min(len(sequence), idx + context_range)
            
            pre_context = sequence[start:idx]
            post_context = sequence[idx+1:end]
            
            if is_valid:
                highlighted_base = f'<span class="mutation-highlight">{alt}</span>'
            else:
                highlighted_base = f'<span class="mutation-highlight" style="background: #ed8936;">{actual_base}</span>'
            
            display_seq = pre_context + highlighted_base + post_context
        
        st.markdown(f"<div class='sequence-container'>{display_seq}</div>", unsafe_allow_html=True)
        
        st.markdown("---")
        
        download_content = create_download_content(sequence, adjusted_position, ref, alt, accession)
        
        col1, col2, col3 = st.columns([1, 1, 1])
        
        with col2:
            st.download_button(
                label="⬇️ Download Sequence",
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
            'effect': analysis.get('effect', 'N/A'),
            'validation': 'PASS' if is_valid else 'FAIL'
        })
    
    if st.session_state.analysis_history:
        st.markdown("---")
        st.subheader("📜 Analysis History")
        df = pd.DataFrame(st.session_state.analysis_history)
        st.dataframe(df, use_container_width=True, hide_index=True)

if __name__ == "__main__":
    main()
