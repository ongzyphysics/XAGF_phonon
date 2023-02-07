function Phonon_out = ReadStoredGFData(fname_in,w_in)
      if ne(exist(fname_in,'file'),2)
          error(sprintf('<!!!> %s does not exist.', fname_read));
      end
      load(fname_in,'Phonon_stored','w_stored');
      if w_in==w_stored
          Phonon_out = Phonon_stored;
      else
          error('<!!!> Stored left lead data is incompatible.');
      end
          
      clear Phonon_stored w_stored;              
end