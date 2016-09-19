#ifndef PATH_CODEC_H
#define PATH_CODEC_H

typedef struct path_encoder_t  {
	unsigned int *buffer;
	int word_offset;
	int bit_offset;
	int bit_width;
} path_encoder_t;

void encoder_init(path_encoder_t &enc, unsigned int *buffer, int bit_width)
{
	enc.buffer = buffer;
	enc.buffer[0] = 0;
	enc.word_offset = 0;
	enc.bit_offset = 0;
	enc.bit_width = bit_width;
}

void encoder_write(path_encoder_t &enc, unsigned int edge_index)
{
	assert((edge_index & ((1 << enc.bit_width)-1)) == edge_index);

	int excess_bit = enc.bit_offset+enc.bit_width - 32;

	if (excess_bit >= 0) {
		//printf("case 1\n");
		enc.buffer[enc.word_offset] |= (edge_index << enc.bit_offset);
		++enc.word_offset;
		enc.buffer[enc.word_offset] = 0 | (edge_index >> (32-enc.bit_offset));
		enc.bit_offset = excess_bit;
	} else {
		//printf("case 2\n");
		enc.buffer[enc.word_offset] |= (edge_index << enc.bit_offset);
		//if (excess_bit == 0) {
			//++enc.word_offset;
			//enc.bit_offset = 0;
			//enc.buffer[enc.word_offset] = 0;
		//} else {
			enc.bit_offset += enc.bit_width;
		//}
	}
}

void encoder_write_trailer(path_encoder_t &enc)
{
	encoder_write(enc, (1 << enc.bit_width)-1);
}

int encoder_get_num_words(path_encoder_t &enc)
{
	return enc.bit_offset == 0 ? enc.word_offset : enc.word_offset+1;
}

typedef struct path_decoder_t  {
	unsigned int *buffer;
	int word_offset;
	int bit_offset;
	int bit_width;
	unsigned int bit_mask;
	bool has_data;
} path_decoder_t;

void decoder_init(path_decoder_t &enc, unsigned int *buffer, int bit_width)
{
	enc.buffer = buffer;
	enc.word_offset = 0;
	enc.bit_offset = 0;
	enc.bit_width = bit_width;
	enc.has_data = true;
	enc.bit_mask = (1 << bit_width) - 1;
}

unsigned int decoder_read(path_decoder_t &enc)
{
	unsigned int val;

	if (!enc.has_data) {
		return enc.bit_mask;
	}

	int excess_bit = enc.bit_offset+enc.bit_width-32;

	if (excess_bit > 0) {
		//printf("case 2\n");

		val = ((enc.buffer[enc.word_offset+1] << (32-enc.bit_offset)) | (enc.buffer[enc.word_offset] >> enc.bit_offset)) & enc.bit_mask;

		++enc.word_offset;
		enc.bit_offset = excess_bit;
	} else {
		//printf("case 3\n");

		val = (enc.buffer[enc.word_offset] >> enc.bit_offset) & enc.bit_mask;

		if (excess_bit == 0) {
			++enc.word_offset;
			enc.bit_offset = 0;
		} else {
			enc.bit_offset += enc.bit_width;
		}
	}

	if (val == enc.bit_mask) {
		enc.has_data = false;
	}

	return val;
}

//bool decoder_has_data(path_decoder_t &enc)
//{
	//return enc.has_data;
//}

#endif
