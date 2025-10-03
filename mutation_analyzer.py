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
    
    .deletion-highlight {
        background: #fc8181;
        color: white;
        padding: 2px 6px;
        border-radius: 4px;
        font-weight: 600;
        text-decoration: line-through;
    }
    
    .insertion-highlight {
        background: #48bb78;
        color: white;
        padding: 2px 6px;
        border-radius: 4px;
        font-weight: 600;
    }
    
    .duplication-highlight {
        background: #4299e1;
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
    
    patterns = {
        'substitution_with_gene': r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):([cg])\.(\d+)([ACGT])>([ACGT])',
        'substitution_with_accession': r'([A-Z0-9_]+(?:\.\d+)?):([cg])\.(\d+)([ACGT])>([ACGT])',
        'intronic_with_gene': r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):c\.(\d+)([\+\-])(\d+)([ACGT])>([ACGT])',
        'intronic_with_accession': r'([A-Z0-9_]+(?:\.\d+)?):c\.(\d+)([\+\-])(\d+)([ACGT])>([ACGT])',
        'deletion_with_gene': r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):([cg])\.(\d+)(?:_(\d+))?del([ACGT]*)',
        'deletion_with_accession': r'([A-Z0-9_]+(?:\.\d+)?):([cg])\.(\d+)(?:_(\d+))?del([ACGT]*)',
        'intronic_deletion_with_gene': r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):c\.(\d+)([\+\-])(\d+)del',
        'intronic_deletion_with_accession': r'([A-Z0-9_]+(?:\.\d+)?):c\.(\d+)([\+\-])(\d+)del',
        'insertion_with_gene': r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):([cg])\.(\d+)_(\d+)ins([ACGT]+)',
        'insertion_with_accession': r'([A-Z0-9_]+(?:\.\d+)?):([cg])\.(\d+)_(\d+)ins([ACGT]+)',
        'duplication_with_gene': r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):([cg])\.(\d+)(?:_(\d+))?dup',
        'duplication_with_accession': r'([A-Z0-9_]+(?:\.\d+)?):([cg])\.(\d+)(?:_(\d+))?dup',
        'delins_with_gene': r'([A-Z0-9_]+(?:\.\d+)?)\(([A-Z0-9]+)\):([cg])\.(\d+)(?:_(\d+))?delins([ACGT]+)',
        'delins_with_accession': r'([A-Z0-9_]+(?:\.\d+)?):([cg])\.(\d+)(?:_(\d+))?delins([ACGT]+)',
    }
    
    match = re.search(patterns['substitution_with_gene'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'coord_type': match.group(3),
            'position': int(match.group(4)),
            'ref': match.group(5).upper(),
            'alt': match.group(6).upper(),
            'variant_type': 'substitution'
        }
    
    match = re.search(patterns['substitution_with_accession'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'coord_type': match.group(2),
            'position': int(match.group(3)),
            'ref': match.group(4).upper(),
            'alt': match.group(5).upper(),
            'variant_type': 'substitution'
        }
    
    match = re.search(patterns['intronic_with_gene'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'coord_type': 'c',
            'position': int(match.group(3)),
            'intronic_sign': match.group(4),
            'intronic_offset': int(match.group(5)),
            'ref': match.group(6).upper(),
            'alt': match.group(7).upper(),
            'variant_type': 'intronic_substitution'
        }
    
    match = re.search(patterns['intronic_with_accession'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'coord_type': 'c',
            'position': int(match.group(2)),
            'intronic_sign': match.group(3),
            'intronic_offset': int(match.group(4)),
            'ref': match.group(5).upper(),
            'alt': match.group(6).upper(),
            'variant_type': 'intronic_substitution'
        }
    
    match = re.search(patterns['deletion_with_gene'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'coord_type': match.group(3),
            'position': int(match.group(4)),
            'end_position': int(match.group(5)) if match.group(5) else int(match.group(4)),
            'deleted_seq': match.group(6).upper() if match.group(6) else None,
            'variant_type': 'deletion'
        }
    
    match = re.search(patterns['deletion_with_accession'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'coord_type': match.group(2),
            'position': int(match.group(3)),
            'end_position': int(match.group(4)) if match.group(4) else int(match.group(3)),
            'deleted_seq': match.group(5).upper() if match.group(5) else None,
            'variant_type': 'deletion'
        }
    
    match = re.search(patterns['intronic_deletion_with_gene'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'coord_type': 'c',
            'position': int(match.group(3)),
            'intronic_sign': match.group(4),
            'intronic_offset': int(match.group(5)),
            'variant_type': 'intronic_deletion'
        }
    
    match = re.search(patterns['intronic_deletion_with_accession'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'coord_type': 'c',
            'position': int(match.group(2)),
            'intronic_sign': match.group(3),
            'intronic_offset': int(match.group(4)),
            'variant_type': 'intronic_deletion'
        }
    
    match = re.search(patterns['insertion_with_gene'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'coord_type': match.group(3),
            'position': int(match.group(4)),
            'end_position': int(match.group(5)),
            'inserted_seq': match.group(6).upper(),
            'variant_type': 'insertion'
        }
    
    match = re.search(patterns['insertion_with_accession'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'coord_type': match.group(2),
            'position': int(match.group(3)),
            'end_position': int(match.group(4)),
            'inserted_seq': match.group(5).upper(),
            'variant_type': 'insertion'
        }
    
    match = re.search(patterns['duplication_with_gene'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'coord_type': match.group(3),
            'position': int(match.group(4)),
            'end_position': int(match.group(5)) if match.group(5) else int(match.group(4)),
            'variant_type': 'duplication'
        }
    
    match = re.search(patterns['duplication_with_accession'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'coord_type': match.group(2),
            'position': int(match.group(3)),
            'end_position': int(match.group(4)) if match.group(4) else int(match.group(3)),
            'variant_type': 'duplication'
        }
    
    match = re.search(patterns['delins_with_gene'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'gene': match.group(2),
            'coord_type': match.group(3),
            'position': int(match.group(4)),
            'end_position': int(match.group(5)) if match.group(5) else int(match.group(4)),
            'inserted_seq': match.group(6).upper(),
            'variant_type': 'delins'
        }
    
    match = re.search(patterns['delins_with_accession'], input_text, re.IGNORECASE)
    if match:
        return {
            'accession': match.group(1),
            'coord_type': match.group(2),
            'position': int(match.group(3)),
            'end_position': int(match.group(4)) if match.group(4) else int(match.group(3)),
            'inserted_seq': match.group(5).upper(),
            'variant_type': 'delins'
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

def extract_cds_info(genbank_record):
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

def cds_to_genomic_position(cds_pos, cds_location, intronic_offset=0):
    if not cds_location:
        return None
    
    parts = cds_location.parts if hasattr(cds_location, 'parts') else [cds_location]
    
    cumulative = 0
    for part in parts:
        part_length = len(part)
        if cumulative + part_length >= cds_pos:
            offset_in_part = cds_pos - cumulative - 1
            genomic_pos = int(part.start) + offset_in_part + intronic_offset
            return genomic_pos
        cumulative += part_length
    
    if parts:
        last_part = parts[-1]
        genomic_pos = int(last_part.end) + intronic_offset
        return genomic_pos
    
    return None

def get_sequence_for_analysis(genbank_record, mutation_data):
    coord_type = mutation_data.get('coord_type', 'c')
    variant_type = mutation_data.get('variant_type')
    position = mutation_data.get('position')
    end_position = mutation_data.get('end_position', position)
    
    cds_seq, cds_location = extract_cds_info(genbank_record)
    
    if coord_type == 'c' and variant_type not in ['intronic_substitution', 'intronic_deletion']:
        if cds_seq:
            return cds_seq, 'CDS', cds_location, position, end_position, None
        else:
            st.warning("CDS not found, using full sequence")
    
    if variant_type in ['intronic_substitution', 'intronic_deletion']:
        intronic_sign = mutation_data.get('intronic_sign')
        intronic_offset = mutation_data.get('intronic_offset', 0)
        
        offset_value = intronic_offset if intronic_sign == '+' else -intronic_offset
        genomic_pos = cds_to_genomic_position(position, cds_location, offset_value)
        
        if genomic_pos:
            position = genomic_pos
            end_position = genomic_pos
            st.info(f"📍 Mapped intronic position to genomic coordinate: {genomic_pos}")
        else:
            st.warning("Could not map intronic position to genomic coordinates")
    
    seq_length = len(genbank_record.seq)
    
    if seq_length > 10000000:
        st.warning(f"⚠️ Large genomic sequence detected ({seq_length:,} bp). Extracting relevant region...")
        
        context = 5000
        start = max(0, position - context)
        end = min(seq_length, max(end_position, position) + context)
        
        try:
            sub_seq = genbank_record.seq[start:end]
            seq_str = get_sequence_safely(sub_seq)
            
            if seq_str:
                adjusted_position = position - start
                adjusted_end = end_position - start
                return seq_str, 'Genomic (Region)', cds_location, adjusted_position, adjusted_end, start
            else:
                st.error("Unable to extract sequence from large genomic record")
                return None, None, None, None, None, None
        except Exception as e:
            st.error(f"Error extracting region: {str(e)}")
            return None, None, None, None, None, None
    
    seq_str = get_sequence_safely(genbank_record.seq)
    if seq_str:
        return seq_str, 'Genomic', cds_location, position, end_position, None
    else:
        st.error("Unable to read sequence from record")
        return None, None, None, None, None, None

def analyze_substitution(sequence, position, ref, alt):
    analysis = {'variant_type': 'Substitution'}
    
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

def analyze_deletion(sequence, position, end_position, expected_seq=None):
    analysis = {'variant_type': 'Deletion'}
    
    start_idx = position - 1
    end_idx = end_position
    
    if start_idx < 0 or end_idx > len(sequence):
        analysis['validation'] = '❌ Position out of range'
        analysis['deleted_length'] = 0
        return analysis
    
    actual_deleted = sequence[start_idx:end_idx]
    analysis['actual_deleted_seq'] = actual_deleted.upper()
    analysis['deleted_length'] = len(actual_deleted)
    
    if expected_seq:
        if actual_deleted.upper() == expected_seq.upper():
            analysis['validation'] = '✅ Deleted sequence matches'
            analysis['validation_status'] = 'match'
        else:
            analysis['validation'] = f'⚠️ Mismatch! Expected {expected_seq}, found {actual_deleted.upper()}'
            analysis['validation_status'] = 'mismatch'
    else:
        analysis['validation'] = '✅ Deletion detected'
        analysis['validation_status'] = 'match'
    
    if analysis['deleted_length'] % 3 == 0:
        analysis['frameshift'] = 'In-frame deletion'
        analysis['effect'] = 'In-frame'
    else:
        analysis['frameshift'] = 'Frameshift deletion'
        analysis['effect'] = 'Frameshift'
    
    context_start = max(0, start_idx - 10)
    context_end = min(len(sequence), end_idx + 10)
    analysis['context'] = sequence[context_start:context_end]
    
    return analysis

def analyze_insertion(sequence, position, end_position, inserted_seq):
    analysis = {'variant_type': 'Insertion'}
    
    start_idx = position - 1
    end_idx = end_position
    
    if start_idx < 0 or end_idx > len(sequence):
        analysis['validation'] = '❌ Position out of range'
        return analysis
    
    flanking_seq = sequence[start_idx:end_idx]
    analysis['flanking_seq'] = flanking_seq.upper()
    analysis['inserted_seq'] = inserted_seq.upper()
    analysis['inserted_length'] = len(inserted_seq)
    analysis['validation'] = '✅ Insertion position valid'
    analysis['validation_status'] = 'match'
    
    if len(inserted_seq) % 3 == 0:
        analysis['frameshift'] = 'In-frame insertion'
        analysis['effect'] = 'In-frame'
    else:
        analysis['frameshift'] = 'Frameshift insertion'
        analysis['effect'] = 'Frameshift'
    
    context_start = max(0, start_idx - 10)
    context_end = min(len(sequence), end_idx + 10)
    analysis['context'] = sequence[context_start:context_end]
    
    return analysis

def analyze_duplication(sequence, position, end_position):
    analysis = {'variant_type': 'Duplication'}
    
    start_idx = position - 1
    end_idx = end_position
    
    if start_idx < 0 or end_idx > len(sequence):
        analysis['validation'] = '❌ Position out of range'
        return analysis
    
    duplicated_seq = sequence[start_idx:end_idx]
    analysis['duplicated_seq'] = duplicated_seq.upper()
    analysis['duplicated_length'] = len(duplicated_seq)
    analysis['validation'] = '✅ Duplication region identified'
    analysis['validation_status'] = 'match'
    
    if len(duplicated_seq) % 3 == 0:
        analysis['frameshift'] = 'In-frame duplication'
        analysis['effect'] = 'In-frame'
    else:
        analysis['frameshift'] = 'Frameshift duplication'
        analysis['effect'] = 'Frameshift'
    
    context_start = max(0, start_idx - 10)
    context_end = min(len(sequence), end_idx + 10)
    analysis['context'] = sequence[context_start:context_end]
    
    return analysis

def analyze_delins(sequence, position, end_position, inserted_seq):
    analysis = {'variant_type': 'Deletion-Insertion'}
    
    start_idx = position - 1
    end_idx = end_position
    
    if start_idx < 0 or end_idx > len(sequence):
        analysis['validation'] = '❌ Position out of range'
        return analysis
    
    deleted_seq = sequence[start_idx:end_idx]
    analysis['deleted_seq'] = deleted_seq.upper()
    analysis['deleted_length'] = len(deleted_seq)
    analysis['inserted_seq'] = inserted_seq.upper()
    analysis['inserted_length'] = len(inserted_seq)
    analysis['validation'] = '✅ Delins region identified'
    analysis['validation_status'] = 'match'
    
    net_change = len(inserted_seq) - len(deleted_seq)
    if net_change % 3 == 0:
        analysis['frameshift'] = 'In-frame delins'
        analysis['effect'] = 'In-frame'
    else:
        analysis['frameshift'] = 'Frameshift delins'
        analysis['effect'] = 'Frameshift'
    
    analysis['net_change'] = net_change
    
    context_start = max(0, start_idx - 10)
    context_end = min(len(sequence), end_idx + 10)
    analysis['context'] = sequence[context_start:context_end]
    
    return analysis

def highlight_variant(sequence, mutation_data, analysis):
    variant_type = mutation_data.get('variant_type')
    position = mutation_data.get('position')
    end_position = mutation_data.get('end_position', position)
    
    idx = position - 1
    end_idx = end_position
    
    if variant_type == 'substitution':
        ref = mutation_data.get('ref')
        alt = mutation_data.get('alt')
        actual_base = analysis.get('actual_base', ref)
        
        if analysis.get('validation_status') == 'match':
            highlighted = sequence[:idx] + f'<span class="mutation-highlight">{alt}</span>' + sequence[idx+1:]
            return highlighted, True
        else:
            highlighted = sequence[:idx] + f'<span class="mutation-highlight" style="background: #ed8936;">{actual_base}</span>' + sequence[idx+1:]
            return highlighted, False
    
    elif variant_type == 'deletion':
        deleted_part = sequence[idx:end_idx]
        highlighted = sequence[:idx] + f'<span class="deletion-highlight">{deleted_part}</span>' + sequence[end_idx:]
        return highlighted, analysis.get('validation_status') == 'match'
    
    elif variant_type == 'insertion':
        inserted_seq = mutation_data.get('inserted_seq')
        highlighted = sequence[:end_idx] + f'<span class="insertion-highlight">[+{inserted_seq}]</span>' + sequence[end_idx:]
        return highlighted, True
    
    elif variant_type == 'duplication':
        duplicated_part = sequence[idx:end_idx]
        highlighted = sequence[:idx] + f'<span class="duplication-highlight">{duplicated_part}</span>' + sequence[idx:end_idx] + f'<span class="duplication-highlight">{duplicated_part}</span>' + sequence[end_idx:]
        return highlighted, True
    
    elif variant_type == 'delins':
        deleted_part = sequence[idx:end_idx]
        inserted_seq = mutation_data.get('inserted_seq')
        highlighted = sequence[:idx] + f'<span class="deletion-highlight">{deleted_part}</span><span class="insertion-highlight">[+{inserted_seq}]</span>' + sequence[end_idx:]
        return highlighted, True
    
    elif variant_type in ['intronic_substitution', 'intronic_deletion']:
        ref = mutation_data.get('ref', '')
        alt = mutation_data.get('alt', '')
        actual_base = analysis.get('actual_base', ref)
        
        if variant_type == 'intronic_substitution':
            highlighted = sequence[:idx] + f'<span class="mutation-highlight">{alt if alt else actual_base}</span>' + sequence[idx+1:]
        else:
            highlighted = sequence[:idx] + f'<span class="deletion-highlight">{actual_base}</span>' + sequence[idx+1:]
        
        return highlighted, True
    
    return sequence, False

def create_download_content(sequence, mutation_data, analysis):
    variant_type = mutation_data.get('variant_type')
    accession = mutation_data.get('accession', 'unknown')
    position = mutation_data.get('position')
    
    if variant_type == 'substitution':
        ref = mutation_data.get('ref')
        alt = mutation_data.get('alt')
        idx = position - 1
        modified_seq = sequence[:idx] + alt + sequence[idx+1:]
        header = f">Mutated_{accession}_pos{position}_{ref}>{alt}"
    
    elif variant_type == 'deletion':
        end_position = mutation_data.get('end_position')
        idx = position - 1
        end_idx = end_position
        modified_seq = sequence[:idx] + sequence[end_idx:]
        header = f">Deleted_{accession}_pos{position}_{end_position}del"
    
    elif variant_type == 'insertion':
        inserted_seq = mutation_data.get('inserted_seq')
        end_idx = mutation_data.get('end_position')
        modified_seq = sequence[:end_idx] + inserted_seq + sequence[end_idx:]
        header = f">Inserted_{accession}_pos{position}_{mutation_data.get('end_position')}ins{inserted_seq}"
    
    elif variant_type == 'duplication':
        end_position = mutation_data.get('end_position')
        idx = position - 1
        end_idx = end_position
        duplicated = sequence[idx:end_idx]
        modified_seq = sequence[:end_idx] + duplicated + sequence[end_idx:]
        header = f">Duplicated_{accession}_pos{position}_{end_position}dup"
    
    elif variant_type == 'delins':
        end_position = mutation_data.get('end_position')
        inserted_seq = mutation_data.get('inserted_seq')
        idx = position - 1
        end_idx = end_position
        modified_seq = sequence[:idx] + inserted_seq + sequence[end_idx:]
        header = f">Delins_{accession}_pos{position}_{end_position}delins{inserted_seq}"
    
    else:
        modified_seq = sequence
        header = f">Modified_{accession}"
    
    content = header + "\n"
    for i in range(0, len(modified_seq), 60):
        content += modified_seq[i:i+60] + "\n"
    
    return content

def main():
    apply_custom_css()
    
    st.markdown("<h1>🧬 Mutation Analyzer</h1>", unsafe_allow_html=True)
    st.markdown("<p class='subtitle'>Comprehensive genomic variant analysis tool</p>", unsafe_allow_html=True)
    
    if 'analysis_history' not in st.session_state:
        st.session_state.analysis_history = []
    
    col1, col2 = st.columns([3, 1])
    
    with col1:
        mutation_input = st.text_input(
            "Enter Variant",
            placeholder="e.g., NM_002440.4(MSH4):c.23C>T or NC_000001.11:g.123456del",
            help="Supports substitutions, deletions, insertions, duplications, and intronic variants"
        )
    
    with col2:
        st.write("")
        st.write("")
        analyze_button = st.button("🔬 Analyze", use_container_width=True)
    
    with st.expander("ℹ️ Supported Variant Types", expanded=False):
        st.markdown("""
        **Substitutions:**
        - `NM_002440.4(MSH4):c.23C>T` - Coding substitution
        - `NC_000001.11:g.123456A>G` - Genomic substitution
        
        **Deletions:**
        - `NM_002440.4:c.100del` - Single base deletion
        - `NM_002440.4:c.100_105del` - Multi-base deletion
        - `NM_002440.4(MSH4):c.2620-10del` - Intronic deletion
        
        **Insertions:**
        - `NM_002440.4:c.100_101insATC` - Insertion between positions
        
        **Duplications:**
        - `NM_002440.4:c.100dup` - Single base duplication
        - `NM_002440.4:c.100_103dup` - Multi-base duplication
        
        **Deletion-Insertions (delins):**
        - `NM_002440.4:c.100_102delinsAGT` - Replace deleted bases with new sequence
        
        **Intronic Variants:**
        - `NM_002440.4(MSH4):c.1906+8G>T` - Intronic substitution
        - `NM_002440.4:c.234-15C>A` - Intronic variant before exon
        """)
    
    if analyze_button and mutation_input:
        with st.spinner("🔍 Parsing variant..."):
            time.sleep(0.2)
            mutation_data = parse_mutation_input(mutation_input)
        
        if not mutation_data:
            st.markdown("<div class='warning-box'>⚠️ Unable to parse variant format. Please check the supported formats above.</div>", unsafe_allow_html=True)
            return
        
        st.markdown("<div class='success-box'>✅ Variant parsed successfully</div>", unsafe_allow_html=True)
        
        accession = mutation_data.get('accession')
        
        with st.spinner(f"📥 Fetching GenBank record for {accession}..."):
            genbank_record = fetch_genbank_record(accession)
        
        if not genbank_record:
            st.error("Failed to fetch GenBank record from NCBI")
            return
        
        st.success(f"✅ Retrieved record: {len(genbank_record.seq):,} bp")
        
        result = get_sequence_for_analysis(genbank_record, mutation_data)
        
        if result[0] is None:
            st.error("Failed to extract sequence")
            return
        
        sequence, seq_type, cds_location, adj_position, adj_end_position, region_start = result
        
        mutation_data['position'] = adj_position
        mutation_data['end_position'] = adj_end_position
        
        if seq_type == 'CDS':
            st.info(f"📍 Working with CDS sequence: {len(sequence):,} bp")
        elif 'Region' in seq_type:
            st.info(f"📍 Extracted genomic region: {len(sequence):,} bp")
        else:
            st.info(f"📍 Working with {seq_type.lower()} sequence: {len(sequence):,} bp")
        
        variant_type = mutation_data.get('variant_type')
        
        if variant_type in ['substitution', 'intronic_substitution']:
            analysis = analyze_substitution(sequence, adj_position, 
                                          mutation_data.get('ref'), 
                                          mutation_data.get('alt'))
        elif variant_type in ['deletion', 'intronic_deletion']:
            analysis = analyze_deletion(sequence, adj_position, adj_end_position, 
                                       mutation_data.get('deleted_seq'))
        elif variant_type == 'insertion':
            analysis = analyze_insertion(sequence, adj_position, adj_end_position, 
                                        mutation_data.get('inserted_seq'))
        elif variant_type == 'duplication':
            analysis = analyze_duplication(sequence, adj_position, adj_end_position)
        elif variant_type == 'delins':
            analysis = analyze_delins(sequence, adj_position, adj_end_position, 
                                     mutation_data.get('inserted_seq'))
        else:
            st.error("Unsupported variant type")
            return
        
        highlighted_seq, is_valid = highlight_variant(sequence, mutation_data, analysis)
        
        validation_status = analysis.get('validation_status', 'unknown')
        if validation_status == 'match':
            st.markdown(f"""<div class='success-box'>
            ✅ <strong>VALIDATED:</strong> {analysis.get('validation', 'Variant validated')}
            </div>""", unsafe_allow_html=True)
        elif validation_status == 'mismatch':
            st.markdown(f"""<div class='validation-highlight'>
            ⚠️ <strong>VALIDATION WARNING:</strong> {analysis.get('validation', 'Mismatch detected')}
            </div>""", unsafe_allow_html=True)
        else:
            st.info(analysis.get('validation', 'Validation status unknown'))
        
        st.markdown("---")
        st.subheader("📊 Variant Analysis")
        
        position = mutation_data.get('position')
        end_position = mutation_data.get('end_position', position)
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            pos_display = f"{position}" if position == end_position else f"{position}-{end_position}"
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Position</div>
                    <div class='metric-value'>{pos_display}</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col2:
            if variant_type == 'substitution':
                change = f"{mutation_data.get('ref')}>{mutation_data.get('alt')}"
            elif variant_type == 'deletion':
                change = f"del{analysis.get('deleted_length', 0)}bp"
            elif variant_type == 'insertion':
                change = f"ins{analysis.get('inserted_length', 0)}bp"
            elif variant_type == 'duplication':
                change = f"dup{analysis.get('duplicated_length', 0)}bp"
            elif variant_type == 'delins':
                change = f"delins"
            else:
                change = variant_type
            
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Change</div>
                    <div class='metric-value'>{change}</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col3:
            mut_type = analysis.get('mutation_type', analysis.get('variant_type', 'N/A'))
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Type</div>
                    <div class='metric-value' style='font-size: 1.2rem;'>{mut_type}</div>
                </div>
            """, unsafe_allow_html=True)
        
        with col4:
            effect = analysis.get('effect', analysis.get('frameshift', 'N/A'))
            st.markdown(f"""
                <div class='metric-card'>
                    <div class='metric-label'>Effect</div>
                    <div class='metric-value' style='font-size: 1.2rem;'>{effect}</div>
                </div>
            """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 🧪 Molecular Details")
            st.markdown("<div class='info-box'>", unsafe_allow_html=True)
            
            if variant_type in ['substitution', 'intronic_substitution']:
                st.markdown(f"<p><strong>Validation:</strong> {analysis['validation']}</p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Actual Base:</strong> {analysis['actual_base']}</p>", unsafe_allow_html=True)
                
                if 'original_codon' in analysis:
                    st.markdown(f"<p><strong>Original Codon:</strong> {analysis['original_codon']}</p>", unsafe_allow_html=True)
                    st.markdown(f"<p><strong>Mutated Codon:</strong> {analysis['mutated_codon']}</p>", unsafe_allow_html=True)
                    st.markdown(f"<p><strong>Codon Position:</strong> {analysis['codon_position']}</p>", unsafe_allow_html=True)
                
                if 'original_aa' in analysis:
                    st.markdown(f"<p><strong>Amino Acid Change:</strong> {analysis['original_aa']} → {analysis['mutated_aa']}</p>", unsafe_allow_html=True)
            
            elif variant_type in ['deletion', 'intronic_deletion']:
                st.markdown(f"<p><strong>Validation:</strong> {analysis['validation']}</p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Deleted Sequence:</strong> <code>{analysis.get('actual_deleted_seq', 'N/A')}</code></p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Deleted Length:</strong> {analysis.get('deleted_length', 0)} bp</p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Frameshift:</strong> {analysis.get('frameshift', 'N/A')}</p>", unsafe_allow_html=True)
            
            elif variant_type == 'insertion':
                st.markdown(f"<p><strong>Inserted Sequence:</strong> <code>{analysis.get('inserted_seq', 'N/A')}</code></p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Inserted Length:</strong> {analysis.get('inserted_length', 0)} bp</p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Frameshift:</strong> {analysis.get('frameshift', 'N/A')}</p>", unsafe_allow_html=True)
            
            elif variant_type == 'duplication':
                st.markdown(f"<p><strong>Duplicated Sequence:</strong> <code>{analysis.get('duplicated_seq', 'N/A')}</code></p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Duplicated Length:</strong> {analysis.get('duplicated_length', 0)} bp</p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Frameshift:</strong> {analysis.get('frameshift', 'N/A')}</p>", unsafe_allow_html=True)
            
            elif variant_type == 'delins':
                st.markdown(f"<p><strong>Deleted Sequence:</strong> <code>{analysis.get('deleted_seq', 'N/A')}</code></p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Inserted Sequence:</strong> <code>{analysis.get('inserted_seq', 'N/A')}</code></p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Net Change:</strong> {analysis.get('net_change', 0)} bp</p>", unsafe_allow_html=True)
                st.markdown(f"<p><strong>Frameshift:</strong> {analysis.get('frameshift', 'N/A')}</p>", unsafe_allow_html=True)
            
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("### 🎯 Sequence Context")
            st.markdown("<div class='info-box'>", unsafe_allow_html=True)
            st.markdown(f"<p><strong>Accession:</strong> {accession}</p>", unsafe_allow_html=True)
            st.markdown(f"<p><strong>Sequence Length:</strong> {len(sequence):,} bp</p>", unsafe_allow_html=True)
            st.markdown(f"<p><strong>Sequence Type:</strong> {seq_type}</p>", unsafe_allow_html=True)
            
            if 'context' in analysis:
                st.markdown(f"<p><strong>Local Context:</strong> <code>...{analysis['context']}...</code></p>", unsafe_allow_html=True)
            
            if cds_location and seq_type == 'CDS':
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
            idx = adj_position - 1
            end_idx = adj_end_position if adj_end_position else adj_position
            
            mid_point = (idx + end_idx) // 2
            start = max(0, mid_point - context_range)
            end = min(len(sequence), mid_point + context_range)
            
            if variant_type == 'substitution':
                pre_context = sequence[start:idx]
                post_context = sequence[idx+1:end]
                if is_valid:
                    highlighted_base = f'<span class="mutation-highlight">{mutation_data.get("alt")}</span>'
                else:
                    highlighted_base = f'<span class="mutation-highlight" style="background: #ed8936;">{analysis.get("actual_base")}</span>'
                display_seq = pre_context + highlighted_base + post_context
            
            elif variant_type == 'deletion':
                pre_context = sequence[start:idx]
                deleted_part = sequence[idx:end_idx]
                post_context = sequence[end_idx:end]
                display_seq = pre_context + f'<span class="deletion-highlight">{deleted_part}</span>' + post_context
            
            elif variant_type == 'insertion':
                pre_context = sequence[start:end_idx]
                post_context = sequence[end_idx:end]
                inserted_seq = mutation_data.get('inserted_seq')
                display_seq = pre_context + f'<span class="insertion-highlight">[+{inserted_seq}]</span>' + post_context
            
            elif variant_type == 'duplication':
                pre_context = sequence[start:idx]
                duplicated_part = sequence[idx:end_idx]
                post_context = sequence[end_idx:end]
                display_seq = pre_context + f'<span class="duplication-highlight">{duplicated_part}</span>' + duplicated_part + f'<span class="duplication-highlight">{duplicated_part}</span>' + post_context
            
            elif variant_type == 'delins':
                pre_context = sequence[start:idx]
                deleted_part = sequence[idx:end_idx]
                post_context = sequence[end_idx:end]
                inserted_seq = mutation_data.get('inserted_seq')
                display_seq = pre_context + f'<span class="deletion-highlight">{deleted_part}</span><span class="insertion-highlight">[+{inserted_seq}]</span>' + post_context
            
            else:
                display_seq = sequence[start:end]
        
        st.markdown(f"<div class='sequence-container'>{display_seq}</div>", unsafe_allow_html=True)
        
        st.markdown("---")
        
        download_content = create_download_content(sequence, mutation_data, analysis)
        
        col1, col2, col3 = st.columns([1, 1, 1])
        
        with col2:
            st.download_button(
                label="⬇️ Download Modified Sequence",
                data=download_content,
                file_name=f"variant_{accession}_{position}.fasta",
                mime="text/plain",
                use_container_width=True
            )
        
        st.session_state.analysis_history.append({
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'input': mutation_input,
            'accession': accession,
            'position': f"{position}" if position == end_position else f"{position}-{end_position}",
            'type': variant_type,
            'effect': analysis.get('effect', analysis.get('frameshift', 'N/A')),
            'validation': 'PASS' if validation_status == 'match' else 'WARN' if validation_status == 'mismatch' else 'N/A'
        })
    
    if st.session_state.analysis_history:
        st.markdown("---")
        st.subheader("📜 Analysis History")
        df = pd.DataFrame(st.session_state.analysis_history)
        st.dataframe(df, use_container_width=True, hide_index=True)

if __name__ == "__main__":
    main()
