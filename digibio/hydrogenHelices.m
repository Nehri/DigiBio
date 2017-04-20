 %protein = pdbread('3bvu.pdb', 'ModelNum', 1);
% 
% 
% 
% %% 
% clc;
% numHelices = length(protein.Helix);
% 
% counts = [0 0 0];
% 
% %for i = 1:numHelices
% for i = 1:2   
%     startLocation = protein.Helix(i).initSeqNum;
%     helixLength = protein.Helix(i).length;
%     endLocation = startLocation + helixLength;
%     
%     for j = startLocation:endLocation
%         
%         donorAcid = find([protein.Model.Atom.resSeq] == j);
%         donorArray = findstr('N', [protein.Model.Atom(donorAcid).AtomName]);
%         
%         donorIndex = donorArray(1);
%         donorAtomNumber = donorAcid(donorIndex);
%         donorAtom = protein.Model.Atom(donorAtomNumber);
%         vecD = [donorAtom.X donorAtom.Y donorAtom.Z];
%         
%         for k = [3:5]
%             if j+k > endLocation
%                 break
%             end
%            
%             acceptorAcid = find([protein.Model.Atom.resSeq] == j+k);
%             acceptorArray = findstr('O', [protein.Model.Atom(acceptorAcid).AtomName]);
%             
%             acceptorIndex = acceptorArray(1);
%             acceptorAtomNumber = acceptorAcid(acceptorIndex);
%             display(acceptorAtomNumber);
%             acceptorAtom = protein.Model.Atom(acceptorAtomNumber);
%             vecA = [acceptorAtom.X acceptorAtom.Y acceptorAtom.Z];
%             
%             
%             antecedentAtomNumber = acceptorAtomNumber - 1;
%             display(antecedentAtomNumber);
%             antecedentAtom = protein.Model.Atom(antecedentAtomNumber);
%             vecB = [antecedentAtom.X antecedentAtom.Y antecedentAtom.Z];
%             
%             hydrogenAcid = find([protein.Model.Atom.resSeq] == j+k);
%             hydrogenArray = findstr('H', [protein.Model.Atom(hydrogenAcid).AtomName]);
%             
%             % Checks back one amino acid more in case no hydrogen was found
%             if isempty(hydrogenArray)
%                 hydrogenAcid = find([protein.Model.Atom.resSeq] == j+k-1);
%                 hydrogenArray = findstr('H', [protein.Model.Atom(hydrogenAcid).AtomName]);
%             end
%             
%             hydrogenIndex = hydrogenArray(1);
%             hydrogenAtomNumber = hydrogenAcid(hydrogenIndex);
%             display(hydrogenAtomNumber);
%             hydrogenAtom = protein.Model.Atom(hydrogenAtomNumber);
%             vecH = [hydrogenAtom.X hydrogenAtom.Y hydrogenAtom.Z];
%             
%             result = hydrogen_analysis(vecD,vecA,vecH,vecB);
%             
%             if result
%                 counts(k-2) = counts(k-2)+1;
%             end
%         end
%     end 
% end
% 
% display(counts(1));
% display(counts(2));
% display(counts(3));


%protein = pdbread('7AHL_B.pdb', 'ModelNum',1);


% %%
% clc;
% numHelices = length(protein.Helix);
%  
% counts = [0 0 0];
% 
% %for i = 1:numHelices
% for i = 1:1
%     startLocation = protein.Helix(i).initSeqNum;
%     helixLength = protein.Helix(i).length;
%     endLocation = startLocation + helixLength-1;
%     
%     for j = startLocation:endLocation
%         donorAcid = find([protein.Model.Atom.resSeq] == j);
%         display(j);
%         donorArray = findstr('N', [protein.Model.Atom(donorAcid).element]);
%         %donorArray = find([protein.Model.Atom(donorAcid).AtomName] == 'N');
%         donorIndex = donorArray(1);
%         donorAtom = donorAcid(donorIndex);
%         donor = protein.Model.Atom(donorAtom);
%         vecD = [donor.X donor.Y donor.Z];
%         
%         for k = 3:5
%             if (j+k) > endLocation 
%                 break;
%             end
%             
%             acceptorAcid = find([protein.Model.Atom.resSeq] == (j+k));
%             acceptorArray = findstr('O', [protein.Model.Atom(acceptorAcid).element]);
%             %acceptorArray = find([protein.Model.Atom(acceptorAcid).AtomName] == 'O');
%             acceptorIndex = acceptorArray(1);
%             acceptorAtom = acceptorAcid(acceptorIndex);
%             acceptor = protein.Model.Atom(acceptorAtom);
%             vecA = [acceptor.X acceptor.Y acceptor.Z];
%             
%             antecedentAtom = acceptorAtom - 1;
%             antecedent = protein.Model.Atom(antecedentAtom);
%             vecB = [antecedent.X antecedent.Y antecedent.Z];
%             
%             %hydrogenAcid = find([protein.Model.Atom.resSeq] == (j+k));
%             hydrogenArray = findstr('H', [protein.Model.Atom(acceptorAcid).element]);
% 
%             % Checks back one amino acid more in case no hydrogen was found
%             %if isempty(hydrogenArray)
%             %    hydrogenAcid = find([protein.Model.Atom.resSeq] == j+k-1);
%             %    hydrogenArray = find([protein.Model.Atom(acceptorAcid).AtomName] == 'H');
%             %end
%             
%             hydrogenIndex = hydrogenArray(1);
%             hydrogenAtom = acceptorAcid(hydrogenIndex);
%             hydro = protein.Model.Atom(hydrogenAtom);
%             vecH = [hydro.X hydro.Y hydro.Z];
%             vecD
%             vecH
%             vecA
%             vecB
%             %result = hydrogen_analysis(vecD,vecA,vecH,vecB);
%             
%             result = hydrogen_analysis(vecD,vecA,vecH,vecB);
%             display(result);
%             if result
%                 display('hit');
%                 counts(k-2) = counts(k-2)+1;
%             end
% 
%         end
%     end
% end

% Helix_a = struct('Name', 'Alpha helix', 'N_T', 14, 'd_ON', [2.99 0.14], 'd_OH', [2.06 0.16], 'Angle_NHO', [155 11], 'Angle_HOC', [147 9], 'Beta', 27, 'Gamma', -18);
% Helix_3_10 = struct('Name', '3_10 helix', 'N_T', 17, 'd_ON', [3.09 0.14], 'd_OH', [2.17 0.16], 'Angle_NHO', [153 10], 'Angle_HOC', [114 10], 'Beta', 55, 'Gamma', -25);
% Sheet_p = struct('Name', 'Parallel beta sheet', 'N_T', 16, 'd_ON', [2.92 0.14], 'd_OH', [1.97 0.15], 'Angle_NHO', [161 9], 'Angle_HOC', [155 11], 'Beta', 0, 'Gamma', -20);
% Sheet_anti_p = struct('Name', 'Anti-parallel beta sheet', 'N_T', 16, 'd_ON', [2.91 0.14], 'd_OH', [1.96 0.16], 'Angle_NHO', [160 10], 'Angle_HOC', [150 12], 'Beta', -10, 'Gamma', 20);
% Structures = [Helix_a, Helix_3_10, Sheet_p, Sheet_anti_p];
% Let D, A, H, B be the donor, the acceptor, the hydrogen, and the
% antecedent.
% 1. |D - A| < 3.5 Ang
% 2. |H - A| < 2.5 Ang
% 3. Angle(DHA) > 90 Deg
% 4. Angle(DAB) > 90 Deg
% 5. Angle(HAB) > 90 Deg
% Criteria = [3.5, 2.5, 90, 90, 90];
% The criteria are ordered according to the ordering of the book.

%% Testing
protein = pdbread('4ug4H.pdb', 'ModelNum', 1);
%%
numHelix = length(protein.Helix);
count = 0;
for i = 1:numHelix % GO THROUGH ALL HELICES
    startPosition = protein.Helix(i).initSeqNum; % Start location within the sequence.
    endPosition = protein.Helix(i).endSeqNum;
    for seqNumber = startPosition:endPosition % Go through each sequence number.
        % Atom Number Range (ANR)
        donorRange = find([protein.Model.Atom.resSeq] == seqNumber);
        
        % Within the ANR, find the Nitrogen (N) atom.
        for atomIndex = donorRange
            
            currentStr = protein.Model.Atom(atomIndex).AtomName;
            
            if strcmp('N', currentStr)
                donorVec = [protein.Model.Atom(atomIndex).X protein.Model.Atom(atomIndex).Y protein.Model.Atom(atomIndex).Z];
            elseif strcmp('H', currentStr)
                    hydrogenVec = [protein.Model.Atom(atomIndex).X protein.Model.Atom(atomIndex).Y protein.Model.Atom(atomIndex).Z];
            end
        end
        
        for j = 3:5
            acceptorRange = find([protein.Model.Atom.resSeq] == (seqNumber + j));
            for acceptorIndex = acceptorRange
                
                currentStr = protein.Model.Atom(acceptorIndex).AtomName;
                
                if strcmp('O', currentStr)
                    acceptorVec = [protein.Model.Atom(acceptorIndex).X protein.Model.Atom(acceptorIndex).Y protein.Model.Atom(acceptorIndex).Z];
                elseif strcmp('C', currentStr)
                    antecedentVec = [protein.Model.Atom(acceptorIndex).X protein.Model.Atom(acceptorIndex).Y protein.Model.Atom(acceptorIndex).Z];
                end
                
            end
            
%             antecedent_idx = acceptor_idx - 1;
            
%             acceptor = [pdbstruct_Model.Model.Atom(acceptor_idx).X pdbstruct_Model.Model.Atom(acceptor_idx).Y pdbstruct_Model.Model.Atom(acceptor_idx).Z];
%             antecedent = [pdbstruct_Model.Model.Atom(antecedent_idx).X pdbstruct_Model.Model.Atom(antecedent_idx).Y pdbstruct_Model.Model.Atom(antecedent_idx).Z];
            
%             hydrogen_r = findstr('H', [pdbstruct_Model.Model.Atom(acceptor_range).AtomName]);
%             hydrogen_idx = acceptor_range(hydrogen_r(1));
%             
%             hydrogen = [pdbstruct_Model.Model.Atom(hydrogen_idx).X pdbstruct_Model.Model.Atom(hydrogen_idx).Y pdbstruct_Model.Model.Atom(hydrogen_idx).Z];
            
            result = hydrogen_analysis(donorVec, acceptorVec, hydrogenVec, antecedentVec);
            
            if result
                display('Got one');
                count = count + 1;
            end
            
        end
    end
display(count);
end