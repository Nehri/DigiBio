
%% Testing
protein = pdbread('4ug4H.pdb', 'ModelNum', 1);
%%
numSheet = length(protein.Sheet);
for i = 1:numSheet % GO THROUGH ALL HELICES
    startPosition = protein.Sheet(i).initSeqNum; % Start location within the sequence.
    endPosition = protein.Sheet(i).endSeqNum;
    orientation = protein.Sheet(i).sense; % 1 if parallel, -1 if anti-parallel, 0 if first strand
    for seqNumber = startPosition:endPosition % Go through each sequence number.
        
        
        % Atom Number Range (ANR)
        donorRange = find([protein.Model.Atom.resSeq] == seqNumber);
        
        % Within the ANR, find the Nitrogen (N) atom.
        for atomIndex = donorRange
            
            currentStr = protein.Model.Atom(atomIndex).AtomName;
            
            if strcmp('N', currentStr)
                donorVec = [protein.Model.Atom(atomIndex).X protein.Model.Atom(atomIndex).Y protein.Model.Atom(atomIndex).Z];
                break;
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
                elseif strcmp('H', currentStr)
                    hydrogenVec = [protein.Model.Atom(acceptorIndex).X protein.Model.Atom(acceptorIndex).Y protein.Model.Atom(acceptorIndex).Z];
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
            end
            
        end
    end
end