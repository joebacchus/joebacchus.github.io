document.body.style.backgroundColor = '';

const glowButton = document.querySelector('#name .glow span');

glowButton.parentElement.addEventListener('mouseenter', () => {
    document.body.style.backgroundColor = '#1a1a1a';
    glowButton.textContent = 'Vanity';
});

glowButton.parentElement.addEventListener('mouseleave', () => {
    document.body.style.backgroundColor = '';
    glowButton.textContent = 'Joe Bacchus George';
});

// Reset background when coming back to page (mobile back button)
window.addEventListener('pageshow', () => {
    document.body.style.backgroundColor = '';
    glowButton.textContent = 'Joe Bacchus George';
});