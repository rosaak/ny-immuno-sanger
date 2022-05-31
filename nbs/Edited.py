from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os
import pandas as pd
from time import perf_counter
from IPython.display import display
import re
import shutil
import statistics
#import PySimpleGUI as sg


abi_sequence_file = "/Users/rp/Desktop/purge/ny-immuno-sanger/from_ny/Sequencing results/VH/"
VH_sequence_folder = "/Users/rp/Desktop/purge/ny-immuno-sanger/from_ny/Templates/VH/"
VL_sequence_folder = "/Users/rp/Desktop/purge/ny-immuno-sanger/from_ny/Templates/VL/"
abspeeq_data = "/Users/rp/Desktop/purge/ny-immuno-sanger/from_ny/14 clones H3.xlsx"
Output_excel_folder_destination = "/Users/rp/Desktop/purge/ny-immuno-sanger/from_ny/Output/"
Output_excel_file_name = "TestJune2022"



def _abi_trim(seq_record):
    start = False  # flag for starting position of trimmed sequence
    segment = 20  # minimum sequence length
    trim_start = 0  # init start index
    cutoff = 0.01  # default cutoff value for calculating base score

    
    
# calculate base score
     score_list = [cutoff - (10 ** (qual / -25.0)) for qual in seq_record.letter_annotations["phred_quality"]]
    # calculate cummulative score
    # if cummulative value < 0, set it to 0
    # first value is set to 0, because of the assumption that
    # the first base will always be trimmed out

    cummul_score = [0]
    for i in range(1, len(score_list)):
        score = cummul_score[-1] + score_list[i]
        if score < 0:
            cummul_score.append(0)
        else:
            cummul_score.append(score)
            if not start:
                # trim_start = value when cummulative score is first > 0
                trim_start = i
                start = True
                # trim_finish = index of highest cummulative score,
                 # marking the end of sequence segment with highest cummulative score
    trim_finish = cummul_score.index(max(cummul_score))
    return seq_record[trim_start:trim_finish]




os.mkdir(Output_excel_folder_destination+'/'+Output_excel_file_name)
column_names = ['Sequence name', 'Matching template', 'Sequence']
matching_sequence_results = pd.DataFrame(columns = column_names)

#make folder for H3 sequences output
os.mkdir(Output_excel_folder_destination+'/'+Output_excel_file_name+'/'+'H3_sorted_sequences')
count = 0
  
#read abspeeq h3 data    
h_data = pd.read_csv(abspeeq_data, sep='\t')

vh_names = []

for template_strand in os.listdir(VH_sequence_folder):
    template_strand = template_strand.split('VH-')[1]
    template_strand = template_strand.split('.gb')[0]
    vh_names.append(template_strand)
    
h3_data = h_data[['name', 'h3_nt']]
h3_names = h3_data['name'].to_list()
h3_nt = h3_data['h3_nt'].to_list()
h3_dict = pd.Series(h3_nt, index=h3_names).to_dict()
h3_matching_sequence_results = pd.DataFrame(columns = ['H3 template sequence', 'Matching well'])

vl_seq_dictionary = {}
template_name_list = []

# Create a new function called alignment

#def alignment:
for filename in os.listdir(abi_sequence_file):  
        if 'ab1' in filename:                 
            if 'VH60' in filename:

                data = ({'Sequence name':filename,'Matching template':'None','Sequence':'Please check manually'})
                matching_sequence_results = matching_sequence_results.append(data, ignore_index=True)
                dna = SeqIO.read(abi_sequence_file+'/'+filename, 'abi')
                
                dna = _abi_trim(dna)
                quality = dna.letter_annotations["phred_quality"]      
                           
                dna = dna.seq
                h3_template_sequence_list = []
                dna_rev = dna.reverse_complement()
                    
                for h3_sequence_key in h3_dict:
                    if str(h3_dict[h3_sequence_key]) in str(dna_rev):
                        for name in vh_names:
                            if name in h3_sequence_key: 
                                try:
                                    os.mkdir(Output_excel_folder_destination+'/'+Output_excel_file_name+
                                             '/'+'H3_sorted_sequences'+'/'+h3_sequence_key)
                                    shutil.copy2(abi_sequence_file+'/'+filename, 
                                                 Output_excel_folder_destination+
                                                 '/'+Output_excel_file_name+'/'+'H3_sorted_sequences'+'/'+h3_sequence_key)

                                except FileExistsError:
                                    shutil.copy2(abi_sequence_file+'/'+filename, 
                                                 Output_excel_folder_destination+'/'
                                                 +Output_excel_file_name+'/'+'H3_sorted_sequences'+'/'+h3_sequence_key)

                                h3_template_sequence_list.append(h3_sequence_key)
                                h3_matching_sequence_results = h3_matching_sequence_results.append(
                                    {'H3 template sequence':h3_sequence_key, 'Matching well':filename}, ignore_index=True)
                    
                for template_strand in os.listdir(VH_sequence_folder):
                    template_name_list.append(template_strand)
                    parsed_template_strand = template_strand.split('VH-')[1]
                    parsed_template_strand = parsed_template_strand.split('.gb')[0]                                                                        
                    for element in h3_template_sequence_list:
                        if parsed_template_strand in element:
                            print(element)
                            template = SeqIO.read(VH_sequence_folder+'/'+template_strand, 'gb')
                            for a in pairwise2.align.globalms(template.seq, dna_rev, 2, -1000, -1000, -1000, 
                                                              penalize_end_gaps = False):
                                    score = int(a[2])
                                    
                                    # quick check score is reasonable
                                    if score >= 700:
                                        
                                        aligned_seq_a = a[0]
                                        aligned_seq_a = aligned_seq_a.replace('-', '')
                                        aligned_seq_a_start = aligned_seq_a[0:6]
                                        aligned_seq_a_length = len(aligned_seq_a) - 6
                                        aligned_seq_a_end = aligned_seq_a[aligned_seq_a_length:len(aligned_seq_a)]
                                        aligned_seq_b = a[1]
                                        
                                        
                                        # check for gaps in template
                                        if '-' not in aligned_seq_b:
                                            #check aligned seq starts with CTCCAC
                                            if 'CTCCAC' in aligned_seq_a_start:
                                                #check aligned seq ends with CTTTCT
                                                if 'CTTTCT' in aligned_seq_a_end:
                                                
                                                    count += 1
                                                    print(filename)
                                                    print(template_strand)
                                                    print('\n')
                                                    print(format_alignment(*a))
                                                    print('\n')
                                                    
                                                    aligned_seq = a[0]
                                                    aligned_seq = aligned_seq.replace('-', '')
                                                    aligned_seq = Seq(aligned_seq)
                                                    aligned_seq_rev = aligned_seq.reverse_complement()
                                                    dna_string = str(dna)
                                                    aligned_seq_rev = str(aligned_seq_rev)
                                                    start = re.search(aligned_seq_rev, dna_string).start()
                                                    end = re.search(aligned_seq_rev, dna_string).end()
                                                    
                                                    quality = quality[start:end]
                                                    quality_score = statistics.mean(quality)
                                                    lowest_quality = min(quality)

                                                    matching_sequence_results.drop(matching_sequence_results.tail(1).index,inplace=True)
                                                    data = ({'Sequence name':filename,'Matching template':parsed_template_strand,'Sequence':'VH60 Ok', 
                                                             'VH mean quality':quality_score, 'VH lowest quality base':lowest_quality})
                                                    matching_sequence_results = matching_sequence_results.append(data, ignore_index=True)
                                                    filename = filename.split('_GATC')[0]
                                                    filename = filename.split('-VH60')[0]
                                                    vl_seq_dictionary[filename] = parsed_template_strand

                                                    print(count)

                                                    break

                                                else:
                                                    pass
                                            else:
                                                pass
                                            
                                        else:
                                            pass
                                        
                                    else:
                                        pass


                                    
for filename in os.listdir(abi_sequence_file):
        if 'ab1' in filename:
            if 'VL' in filename:
                parsed_filename = filename.split('-VL79')[0]
                if parsed_filename in vl_seq_dictionary.keys():
                    
                    matching_template_strand = vl_seq_dictionary[parsed_filename]
                    dna = SeqIO.read(abi_sequence_file+'/'+filename, 'abi')
                    dna = _abi_trim(dna)
                    
                    quality = dna.letter_annotations["phred_quality"]
                               
                    dna = dna.seq            
                    dna_for = dna.reverse_complement()
                        
                    for template_strand in os.listdir(VL_sequence_folder):
                        if matching_template_strand in template_strand:
                            template = SeqIO.read(VL_sequence_folder+'/'+template_strand, 'gb')
                            for a in pairwise2.align.globalms(template.seq, dna_for, 2, -1000, -1000, -1000, penalize_end_gaps = False):
                                score = int(a[2])
                                
                                #quick check score is reasonable
                                if score >= 700:
                                    
                                    aligned_seq_a = a[0]
                                    aligned_seq_a = aligned_seq_a.replace('-', '')
                                    aligned_seq_a_start = aligned_seq_a[0:6]
                                    aligned_seq_a_length = len(aligned_seq_a) - 6
                                    aligned_seq_a_end = aligned_seq_a[aligned_seq_a_length:len(aligned_seq_a)]
                                    aligned_seq_b = a[1]
                                    
                                    #check for gaps in template during alignment
                                    if '-' not in aligned_seq_b:
                                        #check aligned seq starts with CTCCAC
                                        if 'CTCCAC' in aligned_seq_a_start:
                                            #check aligned seq ends with GCTTGG
                                            if 'GCTTGG' in aligned_seq_a_end:
                                                
                                                
                                                
                                                count += 1
                                                print(filename)
                                                print(template_strand)
                                                print('\n')
                                                print(format_alignment(*a))

                                                print('\n')

                                                aligned_seq = a[0]
                                                aligned_seq = aligned_seq.replace('-', '')
                                                aligned_seq = Seq(aligned_seq)
                                                aligned_seq_rev = aligned_seq.reverse_complement()
                                                dna_string = str(dna)
                                                aligned_seq_rev = str(aligned_seq_rev)
                                                start = re.search(aligned_seq_rev, dna_string).start()
                                                end = re.search(aligned_seq_rev, dna_string).end()
                                                    
                                                quality = quality[start:end]
                                                quality_score = statistics.mean(quality)
                                                lowest_quality = min(quality)
                                                
                                                



                                                data = ({'Sequence name':filename,'Matching template':template_strand,'Sequence':'VL79 Ok', 
                                                         'VL mean quality':quality_score, 'VL lowest quality base':lowest_quality})
                                                matching_sequence_results = matching_sequence_results.append(data, ignore_index=True)

                                                print(count)

                                                break
                                                
                                            else:
                                                pass
                                        else:
                                            pass
                                    else: 
                                        pass
                                else:
                                    pass
                else:     
                    pass
                matching_sequence_results = matching_sequence_results.set_index('Sequence name')





vl_79_ok_list_dictionary = {}
vl_79_mean_quality_dictionary = {}
vl_79_lowest_quality_dictionary = {}


for i in matching_sequence_results.index:
    if 'VL79' in i:
        if 'VL79 Ok' in matching_sequence_results.at[i, 'Sequence']:
            parsed_i = i.split('-VL79')[0]
            vl_mean_quality = matching_sequence_results.at[i, 'VL mean quality']
            vl_lowest_quality = matching_sequence_results.at[i, 'VL lowest quality base']
            matching_template = matching_sequence_results.at[i, 'Matching template']
            vl_79_ok_list_dictionary[parsed_i] = matching_template
            vl_79_mean_quality_dictionary[parsed_i] = vl_mean_quality
            vl_79_lowest_quality_dictionary [parsed_i] = vl_lowest_quality
            
parsed_output = pd.DataFrame(columns=['Matching template', 'Sequenced well', 'VH60', 'VL79', 'VH mean quality', 
                                      'VH lowest quality base', 'VL mean quality', 'VL lowest quality base'])

for i in matching_sequence_results.index:
    if 'VH60' in i:
        if 'VH60 Ok' in matching_sequence_results.at[i, 'Sequence']:
            match_temp_name = matching_sequence_results.at[i, 'Matching template']
            vh_mean_qual = matching_sequence_results.at[i, 'VH mean quality']
            vh_lowest_qual = matching_sequence_results.at[i, 'VH lowest quality base']
            for element in vl_79_ok_list_dictionary.keys():
                if element in i:
                    if match_temp_name in vl_79_ok_list_dictionary[element]:
                        vl_mean_qual = vl_79_mean_quality_dictionary[element]
                        vl_lowest_qual = vl_79_lowest_quality_dictionary[element]
                        mean_error_probability = (10**(-vh_lowest_qual/10) + 10**(-vl_lowest_qual/10))/2
                        parsed_output = parsed_output.append({'Matching template':match_temp_name, 'Sequenced well':i, 'VH60':'Ok', 
                                                              'VL79':'Ok', 'VH mean quality':vh_mean_qual, 'VH lowest quality base':vh_lowest_qual,
                                                              'VL mean quality':vl_mean_qual, 'VL lowest quality base':vl_lowest_qual, 
                                                              'Mean probability of error for lowest quality bases in VH and VL':mean_error_probability}, 
                                                             ignore_index=True)  


                        
parsed_output = parsed_output.sort_values('Matching template', ascending=False)

parsed_output = parsed_output.reset_index(drop=True)

final_output = pd.DataFrame(columns=['Matching template', 'Sequenced well', 'VH60', 'VL79', 'VH lowest quality base', 'VL lowest quality base'])

for i in parsed_output.index:
    current_matching_temp = parsed_output.at[i, 'Matching template']
    seq_name = parsed_output.at[i, 'Sequenced well']
    vh_qual = parsed_output.at[i, 'VH lowest quality base']
    vl_qual = parsed_output.at[i, 'VL lowest quality base']
    mean_error_probability = parsed_output.at[i, 'Mean probability of error for lowest quality bases in VH and VL']
    
    if i >= 1:
        prev_matching_temp = parsed_output.at[i-1, 'Matching template']
        if prev_matching_temp in current_matching_temp:
            final_output = final_output.append({'Matching template':'', 'Sequenced well':seq_name, 'VH60':'Ok', 'VL79':'Ok', 
                                                'VH lowest quality base':vh_qual, 'VL lowest quality base':vl_qual, 
                                                'Mean probability of error for lowest quality bases in VH and VL':mean_error_probability}, ignore_index=True) 
        else:
            final_output = final_output.append({'Matching template':'', 'Sequenced well':'', 'VH60':'', 'VL79':'',
                                                'VH lowest quality base':'', 'VL lowest quality base':'', 
                                               'Mean probability of error for lowest quality bases in VH and VL':''}, ignore_index=True)
            final_output = final_output.append({'Matching template':current_matching_temp, 'Sequenced well':seq_name, 
                                                'VH60':'Ok', 'VL79':'Ok', 'VH lowest quality base':vh_qual, 'VL lowest quality base':vl_qual, 
                                               'Mean probability of error for lowest quality bases in VH and VL':mean_error_probability}, ignore_index=True)

    else:
        final_output = final_output.append({'Matching template':current_matching_temp, 'Sequenced well':seq_name, 
                                            'VH60':'Ok', 'VL79':'Ok', 'VH lowest quality base':vh_qual, 'VL lowest quality base':vl_qual, 
                                           'Mean probability of error for lowest quality bases in VH and VL':mean_error_probability}, ignore_index=True)
        
parsed_output = parsed_output.sort_values('Mean probability of error for lowest quality bases in VH and VL', ascending=True)

parsed_output = parsed_output.reset_index(drop=True)

best_sequences_output = pd.DataFrame(columns=['Matching template', 'Sequenced well', 'VH60', 'VL79', 
                                              'VH lowest quality base', 'VL lowest quality base', 
                                             'Mean probability of error for lowest quality bases in VH and VL'])

for i in parsed_output.index:
    current_matching_temp = parsed_output.at[i, 'Matching template']
    seq_name = parsed_output.at[i, 'Sequenced well']
    vh_qual = parsed_output.at[i, 'VH lowest quality base']
    vl_qual = parsed_output.at[i, 'VL lowest quality base']
    mean_error_probability = parsed_output.at[i, 'Mean probability of error for lowest quality bases in VH and VL']
    
    
    
    if i >= 1:
        existing_sequences = best_sequences_output['Matching template'].to_list()
        if not current_matching_temp in existing_sequences:
            best_sequences_output = best_sequences_output.append({'Matching template':'', 'Sequenced well':'', 'VH60':'', 'VL79':'',
                                                                  'VH lowest quality base':'', 'VL lowest quality base':'', 
                                                                 'Mean probability of error for lowest quality bases in VH and VL':''}, ignore_index=True)
            best_sequences_output = best_sequences_output.append({'Matching template':current_matching_temp, 
                                                                  'Sequenced well':seq_name, 'VH60':'Ok', 'VL79':'Ok', 'VH lowest quality base':vh_qual, 
                                                                  'VL lowest quality base':vl_qual, 
                                                                 'Mean probability of error for lowest quality bases in VH and VL':mean_error_probability},
                                                                 ignore_index=True)
    
    else:
        best_sequences_output = best_sequences_output.append({'Matching template':current_matching_temp, 'Sequenced well':seq_name, 
                                                              'VH60':'Ok', 'VL79':'Ok', 'VH lowest quality base':vh_qual,
                                                              'VL lowest quality base':vl_qual, 
                                                              'Mean probability of error for lowest quality bases in VH and VL':mean_error_probability}, ignore_index=True)
seq_found = best_sequences_output
        
seq_found_list = seq_found['Matching template'].to_list()

for i in seq_found_list:
    try:
        seq_found_list.remove('')

    except ValueError:
        break
        
seq_not_found = []

for i in template_name_list:
    found = False
       
    for element in seq_found_list:
        if element in i:
            found = True
    
    if found == False:
        seq_not_found.append(i)
        seq_found_list.append(i)

seq_not_found_df = pd.DataFrame()

seq_not_found_df['Not found'] = seq_not_found
    
h3_parsed_output = pd.DataFrame()
            
for i in h3_matching_sequence_results.index:
    current_matching_temp = h3_matching_sequence_results.at[i, 'H3 template sequence']
    seq_name = h3_matching_sequence_results.at[i, 'Matching well']

    if i >= 1:
        prev_matching_temp = h3_matching_sequence_results.at[i-1, 'H3 template sequence']
        if prev_matching_temp in current_matching_temp:
            h3_parsed_output = h3_parsed_output.append({'H3 template sequence':'', 'Matching well':seq_name}, ignore_index=True)
        else:
            h3_parsed_output = h3_parsed_output.append({'H3 template sequence':'', 'Matching well':''}, ignore_index=True)
            h3_parsed_output = h3_parsed_output.append({'H3 template sequence':current_matching_temp, 'Matching well':seq_name}, ignore_index=True)

    else:
        h3_parsed_output = h3_parsed_output.append({'H3 template sequence':current_matching_temp, 'Matching well':seq_name}, ignore_index=True)
          
            
final_output.to_excel(Output_excel_folder_destination+'/'+Output_excel_file_name+'/'+'full_sequence_match_'+
                      Output_excel_file_name, engine='xlsxwriter')
best_sequences_output.to_excel(Output_excel_folder_destination+'/'+Output_excel_file_name
                               +'/'+'best_matching_sequences'+Output_excel_file_name, engine='xlsxwriter')
h3_parsed_output.to_excel(Output_excel_folder_destination+'/'+Output_excel_file_name
                          +'/'+'h3_match_'+Output_excel_file_name, engine='xlsxwriter')
seq_not_found_df.to_excel(Output_excel_folder_destination+'/'+Output_excel_file_name
                          +'/'+'sequences_not_found_'+Output_excel_file_name, engine='xlsxwriter')